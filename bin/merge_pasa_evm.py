#!/usr/bin/env python3
"""Merge PASA-updated gene models with EVM genes that PASA didn't touch.

PASA annotation comparison only outputs genes with transcript evidence overlap.
EVM genes without transcript support are silently dropped. This script merges
the PASA output with the original EVM genes, keeping PASA-updated versions
where available and adding back EVM genes that PASA missed.

When PASA has structurally altered an EVM gene (chimeric merge or significant
CDS change), the conflict is arbitrated in two stages:

1. **Exact RNA junction chains** (when ``--splice`` is given): each side's
   primary CDS intron chain is scored against the trusted splice set. A side
   wins outright when its supported-intron fraction beats the other by at
   least ``--junction-delta`` (default 0.25). Splice coordinates are the one
   evidence that speaks directly to intron-chain accuracy; a protein hit can
   be excellent while every donor/acceptor is wrong. Measured on the 2026-09
   rebuild fleet, the DIAMOND-only policy left ~80% of conflicts as ties
   (Osa: 7,212 of 8,979) and every tie kept PASA, which is where the
   EVM->PREFILTER intron-chain precision loss (-10.6 to -13.3 pt) came from.
2. **DIAMOND vs SwissProt** for junction-indeterminate conflicts (both
   single-exon, no splice file, or support fractions within the delta).
   Best hit is lowest E-value, alignment length as tie-breaker.

- If PASA wins or ties: keep PASA gene (preserves isoforms/UTRs)
- If EVM wins: use EVM as the gene, graft PASA isoforms that fall within
  EVM gene boundaries as additional mRNA isoforms. When a splice set is
  loaded, a multi-exon PASA isoform is grafted only if every intron it adds
  is junction-supported — an unsupported grafted chain is a precision hit at
  intron-chain level, the exact currency this arbitration is buying back.
- If no overlap with PASA: add EVM gene as-is (rescued)

Usage:
    merge_pasa_evm.py pasa.gff3 evm.gff3 output.gff3 [genome.fa swissprot.fa]
        [--splice trusted_splice] [--junction-delta 0.25]
"""
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

from refine_boundaries import get_introns, load_splice_junctions


def parse_genes(gff3_path):
    """Parse GFF3, return list of (chrom, start, end, strand, gene_id, [lines])."""
    genes = []
    current_gene = None
    current_lines = []

    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            if parts[2] == "gene":
                if current_gene:
                    genes.append(current_gene + (current_lines,))
                gene_id = re.search(r"(?:^|;)ID=([^;]+)", parts[8])
                gid = gene_id.group(1) if gene_id else f"unknown_{len(genes)}"
                current_gene = (parts[0], int(parts[3]), int(parts[4]), parts[6], gid)
                current_lines = [line]
            elif current_gene:
                current_lines.append(line)

    if current_gene:
        genes.append(current_gene + (current_lines,))

    return genes


def primary_cds_intervals(lines):
    """Return the first transcript's CDS intervals as a coordinate tuple.

    PASA can move splice sites without changing the number of CDS rows.  CDS
    counts therefore cannot distinguish a UTR-only update from replacement of
    the EVM intron chain.  Keep the first-isoform convention used downstream by
    refine_boundaries.py, and compare the actual structure instead.
    """
    first_mrna = None
    intervals = []
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        if parts[2] in ("mRNA", "transcript"):
            match = re.search(r"(?:^|;)ID=([^;\s]+)", parts[8])
            if first_mrna is None and match:
                first_mrna = match.group(1)
        elif parts[2] == "CDS":
            parent = re.search(r"(?:^|;)Parent=([^;\s]+)", parts[8])
            parents = parent.group(1).split(",") if parent else []
            if first_mrna and parents and first_mrna not in parents:
                continue
            intervals.append((int(parts[3]), int(parts[4])))
    return tuple(sorted(set(intervals)))


def get_mrna_ids(lines):
    """Get mRNA IDs from gene's GFF3 lines."""
    mrna_ids = []
    for line in lines:
        parts = line.strip().split("\t")
        if len(parts) >= 9 and parts[2] == "mRNA":
            m = re.search(r"(?:^|;)ID=([^;]+)", parts[8])
            if m:
                mrna_ids.append(m.group(1))
    return mrna_ids


def get_mrna_spans(lines):
    """Get mRNA (start, end, lines) from gene's GFF3 lines.

    Returns list of (mrna_start, mrna_end, mrna_id, [mrna_line + child_lines]).
    """
    mrnas = []
    current_mrna = None
    current_mrna_lines = []

    for line in lines:
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        if parts[2] == "gene":
            continue
        if parts[2] == "mRNA":
            if current_mrna:
                mrnas.append(current_mrna + (current_mrna_lines,))
            mrna_id = re.search(r"(?:^|;)ID=([^;]+)", parts[8])
            mid = mrna_id.group(1) if mrna_id else "unknown"
            current_mrna = (int(parts[3]), int(parts[4]), mid)
            current_mrna_lines = [line]
        elif current_mrna:
            current_mrna_lines.append(line)

    if current_mrna:
        mrnas.append(current_mrna + (current_mrna_lines,))

    return mrnas


def block_ids(lines):
    """Return every ID defined by rows in one gene block."""
    ids = set()
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        match = re.search(r"\bID=([^;\s]+)", fields[8])
        if match:
            ids.add(match.group(1))
    return ids


def disambiguate_block(lines, taken):
    """Rename a rescued block's IDs when they collide with already-emitted ones.

    The rescue test asks whether an EVM gene overlaps a PASA gene *in the same
    (chrom, strand) bucket* (issue #51). PASA reuses ``evm.*`` identifiers, so a
    model PASA moved to the other strand -- or more than 50 kb away -- is
    invisible to that test: the EVM gene looks uncovered, gets rescued, and the
    same ID is written twice at two loci. That is invalid GFF3, and the AGAT
    pass at the end of the PREFILTER rule resolves the collision by *fusing*
    the two records into one gene whose exons sit on both strands. Osa carried
    2,522 such collisions -- 1,261 of them strand-flipped -- and the resulting
    chimeras made EVM's gff3_file_to_proteins.pl emit sequences longer than the
    chromosome (issue #50), stalling the filter phase.

    Renaming keeps both structures and makes the file valid; dropping either
    one would discard a real gene model.

    Returns (lines, ids_now_in_use).
    """
    local = block_ids(lines)
    if not local & taken:
        return lines, local

    suffix = ".rescued"
    n = 1
    while any(f"{i}{suffix}" in taken for i in local):
        n += 1
        suffix = f".rescued{n}"
    mapping = {i: f"{i}{suffix}" for i in local}

    def _rename_id(match):
        return mapping.get(match.group(1), match.group(1))

    def _rename_parents(match):
        return ",".join(mapping.get(p, p) for p in match.group(1).split(","))

    out = []
    for line in lines:
        if line.startswith("#"):
            out.append(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            out.append(line)
            continue
        attrs = re.sub(r"(?<=\bID=)([^;\s]+)", _rename_id, fields[8])
        attrs = re.sub(r"(?<=\bParent=)([^;\s]+)", _rename_parents, attrs)
        fields[8] = attrs
        out.append("\t".join(fields) + "\n")
    return out, set(mapping.values())


def overlap_fraction(s1, e1, s2, e2):
    """Fraction of gene1 covered by gene2."""
    overlap = max(0, min(e1, e2) - max(s1, s2) + 1)
    length = e1 - s1 + 1
    return overlap / length if length > 0 else 0


def batch_extract_proteins(gene_groups, genome_fa, tmpdir, prefix):
    """Extract proteins for multiple genes using gffread.

    gene_groups: dict of {label: gff_lines}
    Returns: dict of {protein_header_id: protein_sequence}
    """
    gff_tmp = os.path.join(tmpdir, f"{prefix}_all.gff3")
    pep_tmp = os.path.join(tmpdir, f"{prefix}_all.pep")

    with open(gff_tmp, "w") as f:
        for lines in gene_groups.values():
            for line in lines:
                f.write(line)

    try:
        result = subprocess.run(
            ["gffread", gff_tmp, "-g", genome_fa, "-y", pep_tmp],
            capture_output=True, timeout=120
        )
    except (subprocess.TimeoutExpired, FileNotFoundError) as exc:
        # Do not swallow the failure silently (issue #17 defect 3): a gffread
        # outage would otherwise make every conflict a no-hit "tie".
        print(f"WARNING: gffread protein extraction failed for '{prefix}': {exc}",
              file=sys.stderr)
        return {}

    if result.returncode != 0:
        err = result.stderr.decode(errors="replace")[:500] if result.stderr else ""
        print(f"WARNING: gffread returned {result.returncode} for '{prefix}': {err}",
              file=sys.stderr)

    if not os.path.exists(pep_tmp) or os.path.getsize(pep_tmp) == 0:
        print(f"WARNING: gffread produced no proteins for '{prefix}' — "
              "conflicts on this side become no-hit.", file=sys.stderr)
        return {}

    proteins = {}
    current_id = None
    current_seq = []
    with open(pep_tmp) as f:
        for line in f:
            if line.startswith(">"):
                if current_id and current_seq:
                    proteins[current_id] = "".join(current_seq).replace("*", "").replace(".", "")
                current_id = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
    if current_id and current_seq:
        proteins[current_id] = "".join(current_seq).replace("*", "").replace(".", "")

    return proteins


def make_diamond_db(protein_fa, tmpdir):
    """Create a DIAMOND database."""
    db_path = os.path.join(tmpdir, "swissprot")
    try:
        subprocess.run(
            ["diamond", "makedb", "--in", protein_fa, "--db", db_path, "--quiet"],
            capture_output=True, timeout=600
        )
        return db_path
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None


def diamond_search(query_fa, db_path, tmpdir, out_name):
    """Run DIAMOND blastp, return dict of {query_id: (evalue, align_len, bitscore)}."""
    out_tmp = os.path.join(tmpdir, f"{out_name}.diamond.tsv")
    try:
        subprocess.run(
            ["diamond", "blastp",
             "--query", query_fa, "--db", db_path, "--out", out_tmp,
             "--outfmt", "6", "qseqid", "evalue", "length", "bitscore",
             "--max-target-seqs", "1", "--evalue", "1e-5",
             "--threads", "2", "--quiet"],
            capture_output=True, timeout=600
        )
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return {}

    results = {}
    if os.path.exists(out_tmp):
        with open(out_tmp) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    qid = parts[0]
                    evalue = float(parts[1])
                    align_len = int(parts[2])
                    bitscore = float(parts[3])
                    if qid not in results or (evalue, -align_len) < (results[qid][0], -results[qid][1]):
                        results[qid] = (evalue, align_len, bitscore)
    return results


def compare_diamond(evm_hit, pasa_hit):
    """Compare two DIAMOND hits. Returns 'evm', 'pasa', or 'tie'.

    Best hit = lowest E-value, alignment length as tie-breaker.
    """
    if evm_hit is None and pasa_hit is None:
        return "tie"
    if evm_hit is None:
        return "pasa"
    if pasa_hit is None:
        return "evm"

    evm_eval, evm_alen, _ = evm_hit
    pasa_eval, pasa_alen, _ = pasa_hit

    # Compare E-value (lower is better) — 10x difference threshold
    if evm_eval < pasa_eval * 0.1:
        return "evm"
    if pasa_eval < evm_eval * 0.1:
        return "pasa"

    # E-values similar — alignment length as tie-breaker (10% threshold)
    if evm_alen > pasa_alen * 1.1:
        return "evm"
    if pasa_alen > evm_alen * 1.1:
        return "pasa"

    return "tie"


def junction_support(cds_intervals, chrom, strand, junctions):
    """Score a primary CDS chain against the trusted splice set.

    Returns (supported, total) over the chain's introns, or None when the
    chain has no intron (single-exon chains carry no intron-chain signal).
    """
    introns = get_introns(sorted(cds_intervals))
    if not introns:
        return None
    pool = junctions.get((chrom, strand), set()) | junctions.get((chrom, "."), set())
    supported = sum(1 for i in introns if i in pool)
    return supported, len(introns)


def arbitrate_by_junctions(evm_cds, pasa_cds, chrom, strand, junctions, delta):
    """Return 'evm', 'pasa', or None (indeterminate) from exact junction support.

    Primary criterion: the side carrying FEWER junction-orphan introns wins.
    A count comparison beats a fraction-delta because one wrong intron in a
    long chain barely moves the fraction (1/6 = 0.17 slips under a 0.25
    delta) while it is exactly what breaks the intron chain. Measured on Osa
    against the reference: of 1,328 conflicts where the EVM chain was correct
    and PASA's was not, 927 (70%) separate by unsupported count and only 3
    more by fraction. Fraction (with ``delta``) stays as the tie-break when
    both sides carry the same number of orphans but different chain lengths.

    Only decides when BOTH chains are multi-exon: a single-exon side has no
    intron chain to score, and pretending absence-of-introns is evidence would
    let RNA-void regions flip on noise. Indeterminate conflicts fall through
    to DIAMOND, preserving the pre-existing behavior for them.
    """
    evm_sup = junction_support(evm_cds, chrom, strand, junctions)
    pasa_sup = junction_support(pasa_cds, chrom, strand, junctions)
    if evm_sup is None or pasa_sup is None:
        return None
    evm_orphans = evm_sup[1] - evm_sup[0]
    pasa_orphans = pasa_sup[1] - pasa_sup[0]
    if pasa_orphans > evm_orphans:
        return "evm"
    if evm_orphans > pasa_orphans:
        return "pasa"
    evm_frac = evm_sup[0] / evm_sup[1]
    pasa_frac = pasa_sup[0] / pasa_sup[1]
    if pasa_frac >= evm_frac + delta:
        return "pasa"
    if evm_frac >= pasa_frac + delta:
        return "evm"
    return None


def isoform_introns_supported(mrna_lines, chrom, strand, junctions):
    """True when every intron of one PASA isoform is in the trusted splice set.

    Single-exon isoforms pass: they add no intron chain, so they cannot hurt
    intron-chain precision.
    """
    intervals = []
    for line in mrna_lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 9 and parts[2] == "CDS":
            intervals.append((int(parts[3]), int(parts[4])))
    introns = get_introns(sorted(set(intervals)))
    if not introns:
        return True
    pool = junctions.get((chrom, strand), set()) | junctions.get((chrom, "."), set())
    return all(i in pool for i in introns)


def build_evm_gene_with_pasa_isoforms(evm_lines, pasa_lines, evm_start, evm_end, evm_gene_id,
                                      chrom=None, strand=None, junctions=None):
    """Build a gene using EVM as base, adding PASA isoforms within EVM boundaries.

    Returns GFF3 lines for the merged gene.
    """
    output_lines = []

    # Write EVM gene and its mRNA/CDS
    for line in evm_lines:
        output_lines.append(line)

    # Find PASA mRNAs that fall within EVM gene boundaries
    pasa_mrnas = get_mrna_spans(pasa_lines)
    grafted = 0
    graft_rejected = 0
    for mrna_start, mrna_end, mrna_id, mrna_lines in pasa_mrnas:
        # Check if this PASA isoform is contained within EVM gene range
        if mrna_start >= evm_start and mrna_end <= evm_end:
            if junctions is not None and not isoform_introns_supported(
                    mrna_lines, chrom, strand, junctions):
                graft_rejected += 1
                continue
            grafted += 1
            for line in mrna_lines:
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                if parts[2] == "mRNA":
                    # Reparent to EVM gene
                    parts[8] = re.sub(r"Parent=[^;]+", f"Parent={evm_gene_id}", parts[8])
                    # Rename to avoid ID collision
                    parts[8] = re.sub(r"ID=([^;]+)", rf"ID=\1_pasa_iso", parts[8])
                    output_lines.append("\t".join(parts) + "\n")
                else:
                    # Reparent children to renamed PASA mRNA
                    old_parent = re.search(r"Parent=([^;]+)", parts[8])
                    if old_parent:
                        parts[8] = re.sub(r"Parent=([^;]+)", rf"Parent=\1_pasa_iso", parts[8])
                    output_lines.append("\t".join(parts) + "\n")

    return output_lines, grafted, graft_rejected


def parse_args(argv):
    """Split optional flags from the positional arguments.

    Flags may appear anywhere; positional order is unchanged from the
    pre-flag interface so existing callers keep working.
    """
    splice_path = None
    junction_delta = 0.25
    positional = []
    i = 0
    while i < len(argv):
        if argv[i] == "--splice" and i + 1 < len(argv):
            splice_path = argv[i + 1]
            i += 2
        elif argv[i] == "--junction-delta" and i + 1 < len(argv):
            junction_delta = float(argv[i + 1])
            i += 2
        else:
            positional.append(argv[i])
            i += 1
    return positional, splice_path, junction_delta


def main():
    positional, splice_path, junction_delta = parse_args(sys.argv[1:])
    if len(positional) < 3:
        print(f"Usage: {sys.argv[0]} pasa.gff3 evm.gff3 output.gff3 [genome.fa swissprot.fa]"
              " [--splice trusted_splice] [--junction-delta 0.25]",
              file=sys.stderr)
        sys.exit(1)

    pasa_gff, evm_gff, out_gff = positional[0], positional[1], positional[2]
    genome_fa = positional[3] if len(positional) > 3 else None
    swissprot_fa = positional[4] if len(positional) > 4 else None

    use_diamond = genome_fa and swissprot_fa

    junctions = None
    if splice_path and os.path.exists(splice_path):
        junctions = load_splice_junctions(splice_path)
        n_junc = sum(len(v) for k, v in junctions.items() if k[1] != ".")
        print(f"Loaded {n_junc} trusted splice junctions for chain arbitration",
              file=sys.stderr)
    elif splice_path:
        print(f"WARNING: splice file '{splice_path}' missing — junction arbitration off,"
              " falling back to DIAMOND-only conflict resolution", file=sys.stderr)

    pasa_genes = parse_genes(pasa_gff)
    evm_genes = parse_genes(evm_gff)

    print(f"Loaded {len(pasa_genes)} PASA genes, {len(evm_genes)} EVM genes", file=sys.stderr)

    # Index PASA genes by (chrom, strand)
    pasa_by_loc = defaultdict(list)
    for i, (chrom, start, end, strand, gid, lines) in enumerate(pasa_genes):
        pasa_by_loc[(chrom, strand)].append((start, end, gid, lines, i))

    # First pass: classify each EVM gene
    rescued_blocks = []  # EVM genes with no PASA coverage, one block per gene
    rescued_count = 0
    conflicts = []  # (evm_idx, evm_info, pasa_idx, pasa_info)

    for ei, (chrom, start, end, strand, gid, lines) in enumerate(evm_genes):
        key = (chrom, strand)
        evm_len = end - start + 1

        best_overlap = 0
        best_pasa = None
        for ps, pe, pgid, plines, pidx in pasa_by_loc.get(key, []):
            frac = overlap_fraction(start, end, ps, pe)
            if frac > best_overlap:
                best_overlap = frac
                best_pasa = (ps, pe, pgid, plines, pidx)

        if best_overlap < 0.5:
            rescued_count += 1
            rescued_blocks.append(lines)
            continue

        # Check for structural alteration
        ps, pe, pgid, plines, pidx = best_pasa
        pasa_len = pe - ps + 1
        is_chimeric = pasa_len > 2 * evm_len
        evm_cds = primary_cds_intervals(lines)
        pasa_cds = primary_cds_intervals(plines)
        # A row-count threshold missed shifted donor/acceptor coordinates and
        # terminal CDS changes whenever PASA retained roughly the same number
        # of segments.  Those models bypassed evidence-based selection and
        # silently replaced the improved EVM chain.  Any coordinate change is
        # structural; DIAMOND still resolves it, with PASA preferred on ties as
        # required by issue #17.
        cds_changed = bool(evm_cds) and evm_cds != pasa_cds

        if is_chimeric or cds_changed:
            conflicts.append((ei, (chrom, start, end, strand, gid, lines),
                              pidx, (ps, pe, pgid, plines), evm_cds, pasa_cds))

    print(f"Found {rescued_count} EVM genes without PASA coverage (rescued)",
          file=sys.stderr)
    print(f"Found {len(conflicts)} structural conflicts to resolve", file=sys.stderr)

    # Second pass: resolve conflicts — exact junction chains first, DIAMOND for
    # the indeterminate remainder.
    # pasa_replaced[pidx] = [evm_block, ...]  — a single PASA gene may be replaced
    # by MULTIPLE EVM genes when a chimeric PASA merge spanned several EVM genes
    # (issue #17 defect 1: a plain dict overwrote all but the last, losing genes).
    pasa_replaced = defaultdict(list)
    stats = {"evm_wins": 0, "pasa_wins": 0, "tie": 0, "grafted_isoforms": 0,
             "junction_evm": 0, "junction_pasa": 0, "graft_rejected": 0}

    def apply_evm_win(evm_info, pidx, pasa_info):
        chrom, start, end, strand, gid, lines = evm_info
        merged_lines, n_grafted, n_rejected = build_evm_gene_with_pasa_isoforms(
            lines, pasa_info[3], start, end, gid,
            chrom=chrom, strand=strand, junctions=junctions)
        pasa_replaced[pidx].append(merged_lines)
        stats["grafted_isoforms"] += n_grafted
        stats["graft_rejected"] += n_rejected

    undecided = conflicts
    if conflicts and junctions is not None:
        undecided = []
        for conflict in conflicts:
            ei, evm_info, pidx, pasa_info, evm_cds, pasa_cds = conflict
            chrom, _, _, strand, _, _ = evm_info
            verdict = arbitrate_by_junctions(evm_cds, pasa_cds, chrom, strand,
                                             junctions, junction_delta)
            if verdict == "evm":
                stats["junction_evm"] += 1
                apply_evm_win(evm_info, pidx, pasa_info)
            elif verdict == "pasa":
                stats["junction_pasa"] += 1     # keep PASA
            else:
                undecided.append(conflict)
        print(f"Junction arbitration: {stats['junction_evm']} EVM, "
              f"{stats['junction_pasa']} PASA, {len(undecided)} indeterminate "
              f"of {len(conflicts)} conflicts (delta={junction_delta})",
              file=sys.stderr)

    if undecided and use_diamond:
        with tempfile.TemporaryDirectory() as tmpdir:
            diamond_db = make_diamond_db(swissprot_fa, tmpdir)
            if not diamond_db:
                print("WARNING: Could not create DIAMOND DB", file=sys.stderr)
                use_diamond = False

            if use_diamond:
                print(f"DIAMOND DB built, comparing {len(undecided)} conflicts...",
                      file=sys.stderr)

                # Batch extract proteins
                evm_gene_map = {}
                pasa_gene_map = {}
                for ci, (ei, evm_info, pidx, pasa_info, _, _) in enumerate(undecided):
                    evm_gene_map[f"evm_{ci}"] = evm_info[5]
                    pasa_gene_map[f"pasa_{ci}"] = pasa_info[3]

                evm_proteins = batch_extract_proteins(evm_gene_map, genome_fa, tmpdir, "evm")
                pasa_proteins = batch_extract_proteins(pasa_gene_map, genome_fa, tmpdir, "pasa")

                # Write query FASTAs
                evm_query = os.path.join(tmpdir, "evm_queries.fa")
                pasa_query = os.path.join(tmpdir, "pasa_queries.fa")

                with open(evm_query, "w") as f:
                    for ci, (_, evm_info, _, _, _, _) in enumerate(undecided):
                        for mid in get_mrna_ids(evm_info[5]):
                            if mid in evm_proteins:
                                f.write(f">evm_{ci}\n{evm_proteins[mid]}\n")
                                break

                with open(pasa_query, "w") as f:
                    for ci, (_, _, _, pasa_info, _, _) in enumerate(undecided):
                        for mid in get_mrna_ids(pasa_info[3]):
                            if mid in pasa_proteins:
                                f.write(f">pasa_{ci}\n{pasa_proteins[mid]}\n")
                                break

                # Run DIAMOND
                evm_hits = diamond_search(evm_query, diamond_db, tmpdir, "evm")
                pasa_hits = diamond_search(pasa_query, diamond_db, tmpdir, "pasa")

                # Resolve each remaining conflict
                for ci, (ei, evm_info, pidx, pasa_info, _, _) in enumerate(undecided):
                    evm_hit = evm_hits.get(f"evm_{ci}")
                    pasa_hit = pasa_hits.get(f"pasa_{ci}")
                    winner = compare_diamond(evm_hit, pasa_hit)

                    if winner == "evm":
                        # EVM clearly better — use EVM gene (preserves CDS
                        # accuracy), graft PASA isoforms within EVM range.
                        stats["evm_wins"] += 1
                        apply_evm_win(evm_info, pidx, pasa_info)
                    else:
                        # PASA wins OR tie -> keep PASA, per the module docstring
                        # (issue #17 defect 2). A "tie" includes the common
                        # both-sides-no-DIAMOND-hit case, so this stops PASA
                        # models being replaced by EVM with zero evidence.
                        if winner == "pasa":
                            stats["pasa_wins"] += 1
                        else:
                            stats["tie"] += 1

    elif undecided and not use_diamond:
        print(f"No DIAMOND DB — keeping PASA for {len(undecided)} unresolved conflicts",
              file=sys.stderr)
        stats["tie"] += len(undecided)

    # Write merged output. Track every ID as it is written so a rescued EVM
    # block that would reuse one is renamed rather than silently duplicated
    # (issue #51 — AGAT fuses same-ID records into strand-mixed chimeras).
    emitted_ids = set()
    renamed_blocks = 0
    with open(out_gff, "w") as f:
        f.write("##gff-version 3\n")
        for i, (_, _, _, _, _, lines) in enumerate(pasa_genes):
            if i in pasa_replaced:
                # One PASA gene may be replaced by several EVM genes (chimeric
                # merge); write every EVM block, not just the last (issue #17).
                for block in pasa_replaced[i]:
                    block, ids = disambiguate_block(block, emitted_ids)
                    emitted_ids |= ids
                    for line in block:
                        f.write(line)
            else:
                # PASA blocks need the same collision check as EVM ones. The
                # check used to run in one direction only -- a PASA block just
                # registered its IDs -- so a PASA gene written late could
                # collide with an EVM block written earlier for a lower
                # coordinate, and nothing renamed it. That is how Ptr_ODB got
                # evm.model.Ptr_Chr01.1056 written twice: PASA carried the ID at
                # 9,947,847-9,948,680 on the minus strand while EVM had it at
                # 9,403,891-9,408,364 on the plus strand, and AGAT fused the two
                # into one 545 kb mRNA with children on both strands.
                lines, ids = disambiguate_block(lines, emitted_ids)
                emitted_ids |= ids
                for line in lines:
                    f.write(line)
        if rescued_blocks:
            f.write("# Rescued EVM genes without PASA coverage\n")
            for block in rescued_blocks:
                deduped, ids = disambiguate_block(block, emitted_ids)
                if deduped is not block:
                    renamed_blocks += 1
                emitted_ids |= ids
                for line in deduped:
                    f.write(line)

    evm_blocks_emitted = sum(len(v) for v in pasa_replaced.values())
    total_genes = len(pasa_genes) - len(pasa_replaced) + evm_blocks_emitted + rescued_count
    print(f"\nResults:", file=sys.stderr)
    print(f"  PASA genes: {len(pasa_genes)} ({len(pasa_replaced)} replaced by "
          f"{evm_blocks_emitted} EVM genes)", file=sys.stderr)
    print(f"  Rescued EVM genes: {rescued_count} "
          f"({renamed_blocks} renamed to avoid an ID collision)", file=sys.stderr)
    print(f"  Total genes: {total_genes}", file=sys.stderr)
    print(f"  Conflicts: {len(conflicts)} — junction EVM: {stats['junction_evm']}, "
          f"junction PASA: {stats['junction_pasa']}, DIAMOND EVM wins: {stats['evm_wins']}, "
          f"DIAMOND PASA wins: {stats['pasa_wins']}, tie: {stats['tie']}", file=sys.stderr)
    print(f"  PASA isoforms grafted into EVM genes: {stats['grafted_isoforms']} "
          f"({stats['graft_rejected']} rejected: unsupported introns)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
