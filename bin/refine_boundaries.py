#!/usr/bin/env python3
"""Refine gene model boundaries using independent evidence.

For each gene in the input GFF3, checks if its CDS boundaries are truncated
by comparing against Helixer, Augustus, and Miniprot models at the same locus.
When strong independent evidence supports a longer/more complete alternative,
the gene model is swapped.

Evidence sources (all independent — no OMArk/BUSCO in selection):
  1. RNA-seq splice junction support (intron validation)
  2. Miniprot protein alignment coverage
  3. Cross-source CDS boundary agreement

Usage:
    python3 refine_boundaries.py input.gff3 \
        --helixer helixer.gff3 \
        --miniprot miniprot.clean.gff3 \
        --splice psiclass.trusted_splice \
        --output refined.gff3 \
        [--augustus augustus.gff3]
"""
import argparse
import re
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# GFF3 parsing
# ---------------------------------------------------------------------------

def parse_genes(gff3_path):
    """Parse GFF3 into gene models with CDS intervals.

    Returns list of dicts with keys:
        chrom, start, end, strand, gene_id, lines, cds_intervals, exon_intervals
    """
    genes = []
    current = None
    current_lines = []

    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            if cols[2] == "gene":
                if current:
                    current["lines"] = current_lines
                    _extract_intervals(current)
                    genes.append(current)
                gid = _get_attr(cols[8], "ID") or f"unknown_{len(genes)}"
                current = {
                    "chrom": cols[0], "start": int(cols[3]), "end": int(cols[4]),
                    "strand": cols[6], "gene_id": gid, "source": cols[1],
                    "lines": [], "cds_intervals": [], "exon_intervals": []
                }
                current_lines = [line]
            elif current:
                current_lines.append(line)

    if current:
        current["lines"] = current_lines
        _extract_intervals(current)
        genes.append(current)

    return genes


def parse_miniprot(gff3_path):
    """Parse miniprot GFF3 (mRNA/CDS only, no gene lines)."""
    models = []
    current = None
    current_lines = []

    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            if cols[2] == "mRNA":
                if current:
                    current["lines"] = current_lines
                    _extract_intervals(current)
                    models.append(current)
                mid = _get_attr(cols[8], "ID") or f"mp_{len(models)}"
                identity = _get_attr(cols[8], "Identity")
                current = {
                    "chrom": cols[0], "start": int(cols[3]), "end": int(cols[4]),
                    "strand": cols[6], "gene_id": mid, "source": "miniprot",
                    "lines": [], "cds_intervals": [], "exon_intervals": [],
                    "identity": float(identity) if identity else 0.0
                }
                # Wrap mRNA in a synthetic gene line
                gene_line = f"{cols[0]}\tminiprot\tgene\t{cols[3]}\t{cols[4]}\t{cols[5]}\t{cols[6]}\t.\tID={mid}_gene\n"
                current_lines = [gene_line, line]
            elif current:
                current_lines.append(line)

    if current:
        current["lines"] = current_lines
        _extract_intervals(current)
        models.append(current)

    return models


def _get_attr(attrs, key):
    m = re.search(rf"{key}=([^;\s]+)", attrs)
    return m.group(1) if m else None


def _extract_intervals(gene):
    """Extract CDS and exon intervals from gene's GFF3 lines.

    Only uses the first mRNA/transcript to avoid mixing isoform structures.
    """
    cds = []
    exons = []
    first_mrna = None
    for line in gene["lines"]:
        cols = line.strip().split("\t")
        if len(cols) < 9:
            continue
        if cols[2] in ("mRNA", "transcript"):
            mid = _get_attr(cols[8], "ID")
            if first_mrna is None:
                first_mrna = mid
            elif mid != first_mrna:
                continue  # skip non-primary isoforms
        elif cols[2] in ("CDS", "exon"):
            parent = _get_attr(cols[8], "Parent")
            if first_mrna and parent and parent != first_mrna:
                continue  # belongs to a different isoform
            if cols[2] == "CDS":
                cds.append((int(cols[3]), int(cols[4])))
            else:
                exons.append((int(cols[3]), int(cols[4])))
    gene["cds_intervals"] = sorted(set(cds))
    gene["exon_intervals"] = sorted(set(exons))


def get_introns(intervals):
    """Derive intron coordinates from sorted exon/CDS intervals."""
    introns = []
    for i in range(len(intervals) - 1):
        intron_start = intervals[i][1] + 1
        intron_end = intervals[i + 1][0] - 1
        if intron_end > intron_start:
            introns.append((intron_start, intron_end))
    return introns


def total_cds_length(cds_intervals):
    return sum(e - s + 1 for s, e in cds_intervals)


# ---------------------------------------------------------------------------
# Splice junction loading
# ---------------------------------------------------------------------------

def load_splice_junctions(splice_path, tolerance=5):
    """Load trusted splice junctions into lookup structure.

    Returns dict[(chrom, strand)] -> list of (start, end) sorted by start.
    Also returns a set for exact matching with tolerance.
    """
    junctions = defaultdict(set)
    if not splice_path:
        return junctions

    with open(splice_path) as f:
        for line in f:
            cols = line.strip().split()  # space or tab separated
            if len(cols) < 5:
                continue
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            strand = cols[4] if cols[4] in ("+", "-") else "."
            # Store exact junction (tolerance handled at query time)
            junctions[(chrom, strand)].add((start, end))
            junctions[(chrom, ".")].add((start, end))

    return junctions


# ---------------------------------------------------------------------------
# Spatial index for overlap queries
# ---------------------------------------------------------------------------

def build_spatial_index(genes):
    """Build dict[(chrom, strand)] -> sorted list of (start, end, index)."""
    idx = defaultdict(list)
    for i, g in enumerate(genes):
        idx[(g["chrom"], g["strand"])].append((g["start"], g["end"], i))
    for key in idx:
        idx[key].sort()
    return idx


def find_overlapping(gene, spatial_idx, models, min_overlap_frac=0.3):
    """Find models overlapping a gene with at least min_overlap_frac reciprocal overlap."""
    key = (gene["chrom"], gene["strand"])
    if key not in spatial_idx:
        return []

    results = []
    gstart, gend = gene["start"], gene["end"]
    gspan = gend - gstart + 1

    for s, e, i in spatial_idx[key]:
        if s > gend:
            break
        if e < gstart:
            continue
        # Calculate overlap
        ov_start = max(gstart, s)
        ov_end = min(gend, e)
        ov_len = ov_end - ov_start + 1
        if ov_len <= 0:
            continue
        mspan = e - s + 1
        ov_frac = min(ov_len / gspan, ov_len / mspan)
        if ov_frac >= min_overlap_frac:
            results.append(models[i])

    return results


# ---------------------------------------------------------------------------
# Scoring functions (independent evidence only)
# ---------------------------------------------------------------------------

def score_splice_support(model, junctions, tolerance=5):
    """Score: fraction of model's introns supported by splice junctions."""
    # Use exon intervals if available, else CDS
    intervals = model["exon_intervals"] if model["exon_intervals"] else model["cds_intervals"]
    introns = get_introns(intervals)

    if not introns:
        return 0.5  # Single-exon: neutral score

    key = (model["chrom"], model["strand"])
    jset = junctions.get(key, set()) | junctions.get((model["chrom"], "."), set())

    supported = 0
    for intron_start, intron_end in introns:
        found = False
        for ds in range(-tolerance, tolerance + 1):
            for de in range(-tolerance, tolerance + 1):
                if (intron_start + ds, intron_end + de) in jset:
                    found = True
                    break
            if found:
                break
        if found:
            supported += 1

    return supported / len(introns)


def score_miniprot_coverage(model, miniprot_idx, miniprot_models):
    """Score: how well miniprot alignments cover this model's CDS."""
    if not model["cds_intervals"]:
        return 0.0

    overlapping = find_overlapping(model, miniprot_idx, miniprot_models, min_overlap_frac=0.2)
    if not overlapping:
        return 0.0

    # Find best miniprot overlap by CDS coverage
    model_cds_bases = set()
    for s, e in model["cds_intervals"]:
        model_cds_bases.update(range(s, e + 1))

    best_cov = 0.0
    for mp in overlapping:
        if mp.get("identity", 0) < 0.4:
            continue
        mp_cds_bases = set()
        for s, e in mp["cds_intervals"]:
            mp_cds_bases.update(range(s, e + 1))
        overlap = len(model_cds_bases & mp_cds_bases)
        cov = overlap / len(model_cds_bases) if model_cds_bases else 0
        best_cov = max(best_cov, cov)

    return best_cov


def score_cross_source(model, other_models, tolerance=10):
    """Score: fraction of CDS boundaries confirmed by another source."""
    if not model["cds_intervals"]:
        return 0.0

    # Collect all boundary coordinates from model
    boundaries = set()
    for s, e in model["cds_intervals"]:
        boundaries.add(("start", s))
        boundaries.add(("end", e))

    if not boundaries:
        return 0.0

    # Collect boundaries from other models
    other_boundaries = set()
    for om in other_models:
        if om["gene_id"] == model["gene_id"]:
            continue
        for s, e in om["cds_intervals"]:
            for t in range(-tolerance, tolerance + 1):
                other_boundaries.add(("start", s + t))
                other_boundaries.add(("end", e + t))

    confirmed = len(boundaries & other_boundaries)
    return confirmed / len(boundaries)


# ---------------------------------------------------------------------------
# Truncation detection
# ---------------------------------------------------------------------------

def is_truncated(sylvan_gene, alternatives):
    """Check if sylvan_gene appears truncated vs alternatives.

    Conservative: requires the alternative to be genuinely longer in CDS.
    """
    if not sylvan_gene["cds_intervals"]:
        return False

    sylvan_cds_len = total_cds_length(sylvan_gene["cds_intervals"])
    sylvan_cds_count = len(sylvan_gene["cds_intervals"])

    for alt in alternatives:
        if not alt["cds_intervals"]:
            continue
        alt_cds_len = total_cds_length(alt["cds_intervals"])

        # Must have more CDS content (at least 20% more)
        if alt_cds_len < sylvan_cds_len * 1.2:
            continue

        # Condition 1: Alternative CDS is >= 1.5x longer
        if alt_cds_len >= sylvan_cds_len * 1.5:
            return True

        # Condition 2: Alternative extends beyond gene boundaries with shared internal structure
        alt_extends = (alt["start"] < sylvan_gene["start"] - 100 or
                       alt["end"] > sylvan_gene["end"] + 100)
        if alt_extends:
            shared = _count_shared_boundaries(sylvan_gene["cds_intervals"],
                                               alt["cds_intervals"], tolerance=10)
            if shared >= 2:
                return True

        # Condition 3: Alternative has >= 2 more CDS AND more CDS content
        if len(alt["cds_intervals"]) >= sylvan_cds_count + 2:
            return True

    return False


def _count_shared_boundaries(cds_a, cds_b, tolerance=10):
    """Count CDS boundaries shared between two models."""
    bounds_a = set()
    for s, e in cds_a:
        bounds_a.add(s)
        bounds_a.add(e)

    count = 0
    for s, e in cds_b:
        for ba in bounds_a:
            if abs(s - ba) <= tolerance or abs(e - ba) <= tolerance:
                count += 1
                break
    return count


# ---------------------------------------------------------------------------
# Gene model replacement
# ---------------------------------------------------------------------------

def rebuild_gene_lines(original, replacement):
    """Rebuild GFF3 lines using replacement's structure but original's gene ID."""
    gene_id = original["gene_id"]
    # Derive mRNA ID from gene_id pattern
    mrna_id = gene_id.replace("evm.TU.", "evm.model.") if "evm.TU." in gene_id else gene_id + ".t1"

    new_lines = []
    chrom = original["chrom"]
    strand = replacement["strand"]

    # Gene line
    gene_start = replacement["start"]
    gene_end = replacement["end"]
    new_lines.append(
        f"{chrom}\tSylvan_refined\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\tID={gene_id}\n"
    )

    # mRNA line
    new_lines.append(
        f"{chrom}\tSylvan_refined\tmRNA\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t"
        f"ID={mrna_id};Parent={gene_id}\n"
    )

    # Exon lines
    exon_intervals = replacement["exon_intervals"] if replacement["exon_intervals"] else replacement["cds_intervals"]
    for i, (s, e) in enumerate(sorted(exon_intervals), 1):
        new_lines.append(
            f"{chrom}\tSylvan_refined\texon\t{s}\t{e}\t.\t{strand}\t.\t"
            f"ID={mrna_id}.exon{i};Parent={mrna_id}\n"
        )

    # CDS lines
    for i, (s, e) in enumerate(sorted(replacement["cds_intervals"]), 1):
        phase = "0"  # Phase would need proper calculation for accuracy
        new_lines.append(
            f"{chrom}\tSylvan_refined\tCDS\t{s}\t{e}\t.\t{strand}\t{phase}\t"
            f"ID=cds.{mrna_id}.{i};Parent={mrna_id}\n"
        )

    return new_lines, gene_start, gene_end


def overlaps_other_gene(candidate, sylvan_genes, original_gene, spatial_idx, overlap_threshold=0.3):
    """Check if candidate model overlaps a DIFFERENT sylvan gene significantly."""
    key = (candidate["chrom"], candidate["strand"])
    if key not in spatial_idx:
        return False

    for s, e, i in spatial_idx[key]:
        if s > candidate["end"]:
            break
        if e < candidate["start"]:
            continue
        other = sylvan_genes[i]
        if other["gene_id"] == original_gene["gene_id"]:
            continue
        ov_start = max(candidate["start"], s)
        ov_end = min(candidate["end"], e)
        ov_len = ov_end - ov_start + 1
        other_span = e - s + 1
        if ov_len > 0 and ov_len / other_span > overlap_threshold:
            return True
    return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Refine gene model boundaries using independent evidence.")
    parser.add_argument("input_gff", help="Input GFF3 (e.g., Sylvan_merged.gff3)")
    parser.add_argument("--helixer", help="Helixer GFF3")
    parser.add_argument("--augustus", help="Augustus GFF3")
    parser.add_argument("--miniprot", help="Miniprot clean GFF3")
    parser.add_argument("--splice", help="Trusted splice junction file (PsiClass format)")
    parser.add_argument("--output", "-o", required=True, help="Output refined GFF3")
    parser.add_argument("--min-score", type=float, default=0.6,
                        help="Minimum composite score to swap (default: 0.6)")
    args = parser.parse_args()

    # Parse input genes
    print(f"Parsing input GFF3: {args.input_gff}", file=sys.stderr)
    sylvan_genes = parse_genes(args.input_gff)
    print(f"  {len(sylvan_genes)} genes loaded", file=sys.stderr)

    # Parse alternative sources
    alt_sources = {}
    if args.helixer:
        print(f"Parsing Helixer: {args.helixer}", file=sys.stderr)
        alt_sources["helixer"] = parse_genes(args.helixer)
        print(f"  {len(alt_sources['helixer'])} genes", file=sys.stderr)
    if args.augustus:
        print(f"Parsing Augustus: {args.augustus}", file=sys.stderr)
        alt_sources["augustus"] = parse_genes(args.augustus)
        print(f"  {len(alt_sources['augustus'])} genes", file=sys.stderr)
    if args.miniprot:
        print(f"Parsing Miniprot: {args.miniprot}", file=sys.stderr)
        alt_sources["miniprot"] = parse_miniprot(args.miniprot)
        print(f"  {len(alt_sources['miniprot'])} models", file=sys.stderr)

    if not alt_sources:
        print("WARNING: No alternative sources provided. Copying input to output.", file=sys.stderr)
        with open(args.input_gff) as fin, open(args.output, "w") as fout:
            fout.write(fin.read())
        return

    # Load splice junctions
    junctions = load_splice_junctions(args.splice)
    n_junctions = sum(len(v) for v in junctions.values())
    print(f"Loaded splice junctions: {n_junctions // 11} unique (with tolerance)", file=sys.stderr)

    # Build spatial indices
    sylvan_idx = build_spatial_index(sylvan_genes)
    alt_indices = {}
    for name, models in alt_sources.items():
        alt_indices[name] = build_spatial_index(models)

    # Build miniprot index separately for coverage scoring
    miniprot_idx = alt_indices.get("miniprot", {})
    miniprot_models = alt_sources.get("miniprot", [])

    # Process each gene
    refined_count = 0
    candidate_count = 0
    output_genes = []  # list of (start, end, lines)

    for gene in sylvan_genes:
        # Collect all overlapping alternative models
        all_alternatives = []
        for name, models in alt_sources.items():
            overlapping = find_overlapping(gene, alt_indices[name], models, min_overlap_frac=0.3)
            for m in overlapping:
                m["_source_name"] = name
            all_alternatives.extend(overlapping)

        # Check if gene appears truncated
        if not is_truncated(gene, all_alternatives):
            output_genes.append((gene["start"], gene["end"], gene["lines"]))
            continue

        candidate_count += 1

        # Score original gene
        orig_splice = score_splice_support(gene, junctions)

        # Score each alternative
        best_alt = None
        best_score = 0.0
        best_scores = None

        for alt in all_alternatives:
            # Skip if it overlaps another Sylvan gene (chimeric)
            if overlaps_other_gene(alt, sylvan_genes, gene, sylvan_idx, overlap_threshold=0.3):
                continue

            # Skip miniprot with low identity
            if alt.get("identity", 1.0) < 0.4:
                continue

            # Must have CDS and be longer than original
            if not alt["cds_intervals"]:
                continue
            if total_cds_length(alt["cds_intervals"]) <= total_cds_length(gene["cds_intervals"]):
                continue

            # Score
            splice_sc = score_splice_support(alt, junctions)
            mp_sc = score_miniprot_coverage(alt, miniprot_idx, miniprot_models) if miniprot_models else 0.0

            # Cross-source: check against other alternatives from different sources
            other_alts = [a for a in all_alternatives
                          if a.get("_source_name") != alt.get("_source_name")]
            cross_sc = score_cross_source(alt, other_alts)

            composite = 0.4 * splice_sc + 0.3 * mp_sc + 0.3 * cross_sc

            # Check minimum thresholds
            scores_above = sum(1 for s in [splice_sc, mp_sc, cross_sc] if s >= 0.3)
            if scores_above < 2:
                continue

            # Don't regress on splice support
            if splice_sc < orig_splice - 0.1:
                continue

            if composite > best_score:
                best_score = composite
                best_alt = alt
                best_scores = (splice_sc, mp_sc, cross_sc)

        if best_alt and best_score >= args.min_score:
            # Only allow Helixer/Augustus as replacement sources (not miniprot)
            # Miniprot gives approximate exon boundaries; RNA-seq splice junctions are authoritative
            if best_alt.get("_source_name") == "miniprot":
                output_genes.append((gene["start"], gene["end"], gene["lines"]))
                continue

            # Swap
            new_lines, new_start, new_end = rebuild_gene_lines(gene, best_alt)
            output_genes.append((new_start, new_end, new_lines))
            refined_count += 1
            print(f"  REFINED {gene['gene_id']}: "
                  f"CDS {total_cds_length(gene['cds_intervals'])}→{total_cds_length(best_alt['cds_intervals'])}bp, "
                  f"source={best_alt.get('_source_name','?')}, "
                  f"scores=splice:{best_scores[0]:.2f}/mp:{best_scores[1]:.2f}/cross:{best_scores[2]:.2f}="
                  f"{best_score:.2f}",
                  file=sys.stderr)
        else:
            output_genes.append((gene["start"], gene["end"], gene["lines"]))

    # Write output sorted by position
    output_genes.sort(key=lambda x: (x[0], x[1]))

    with open(args.output, "w") as fout:
        fout.write("##gff-version 3\n")
        fout.write(f"# Refined by refine_boundaries.py: {refined_count} genes updated "
                   f"out of {candidate_count} candidates ({len(sylvan_genes)} total)\n")
        for _, _, lines in output_genes:
            for line in lines:
                fout.write(line if line.endswith("\n") else line + "\n")

    print(f"\n=== Refinement Summary ===", file=sys.stderr)
    print(f"Total genes: {len(sylvan_genes)}", file=sys.stderr)
    print(f"Candidates (truncated): {candidate_count}", file=sys.stderr)
    print(f"Refined: {refined_count}", file=sys.stderr)
    print(f"Output: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
