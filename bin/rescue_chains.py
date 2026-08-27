#!/usr/bin/env python3
"""Replace EVM loci with RNA-better single-track gene models.

The primary signal is the fraction of a transcript's CDS introns found at the
exact same coordinates (tolerance zero) in reproducible, canonical STAR
SJ.out.tab records.  Replacement is deliberately locus-atomic and guarded by
complete-ORF, one-locus-only, and exon-growth checks.
"""
import argparse
import os
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field


def attrs(text):
    result = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.strip().split("=", 1)
            result[key] = value
    return result


@dataclass
class Feature:
    raw: str
    seqid: str
    source: str
    kind: str
    start: int
    end: int
    strand: str
    phase: str
    attr: dict


@dataclass
class Gene:
    gene_id: str
    rows: list = field(default_factory=list)
    transcripts: dict = field(default_factory=lambda: defaultdict(list))
    seqid: str = ""
    strand: str = ""
    start: int = 0
    end: int = 0


def parse_feature(line, path, lineno):
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 9:
        raise ValueError(f"{path}:{lineno}: expected 9 GFF3 columns")
    try:
        start, end = int(fields[3]), int(fields[4])
    except ValueError as exc:
        raise ValueError(f"{path}:{lineno}: invalid coordinates") from exc
    return Feature(line.rstrip("\n"), fields[0], fields[1], fields[2],
                   min(start, end), max(start, end), fields[6], fields[7],
                   attrs(fields[8]))


def load_gff(path):
    """Load complete gene bundles; reject ambiguous parentage."""
    headers, features = [], []
    with open(path) as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                headers.append(line.rstrip("\n"))
                continue
            features.append(parse_feature(line, path, lineno))

    genes, tx_gene = {}, {}
    for feat in features:
        if feat.kind == "gene":
            gid = feat.attr.get("ID")
            if not gid or gid in genes:
                raise ValueError(f"{path}: gene has missing/duplicate ID {gid!r}")
            genes[gid] = Gene(gid, seqid=feat.seqid, strand=feat.strand,
                              start=feat.start, end=feat.end)
    for feat in features:
        if feat.kind in ("mRNA", "transcript"):
            tid, parents = feat.attr.get("ID"), feat.attr.get("Parent", "").split(",")
            if not tid or len(parents) != 1 or parents[0] not in genes:
                raise ValueError(f"{path}: transcript {tid!r} lacks one known gene parent")
            tx_gene[tid] = parents[0]
    for feat in features:
        if feat.kind == "gene":
            gid = feat.attr["ID"]
        elif feat.kind in ("mRNA", "transcript"):
            gid = tx_gene[feat.attr["ID"]]
        else:
            parents = feat.attr.get("Parent", "").split(",")
            gids = {tx_gene[p] for p in parents if p in tx_gene}
            if len(gids) != 1:
                raise ValueError(f"{path}: child has ambiguous/unknown Parent: {feat.raw}")
            gid = gids.pop()
            for parent in parents:
                if parent in tx_gene:
                    genes[gid].transcripts[parent].append(feat)
        genes[gid].rows.append(feat)
    return headers, list(genes.values())


def cds_introns(rows):
    cds = sorted((f.start, f.end) for f in rows if f.kind == "CDS")
    return tuple((left[1] + 1, right[0] - 1)
                 for left, right in zip(cds, cds[1:])
                 if left[1] + 1 <= right[0] - 1)


def complete_orf(rows):
    """Require explicit completeness or a coherent full CDS representation."""
    tx = next((f for f in rows if f.kind in ("mRNA", "transcript")), None)
    combined = tx.attr if tx else {}
    bad = (combined.get("valid_ORF", "").lower() == "false" or
           any(combined.get(k, "").lower() == "true" for k in
               ("missing_start_codon", "missing_stop_codon",
                "inframe_stop_codon", "partial_mapping")))
    if bad:
        return False
    if combined.get("valid_ORF", "").lower() == "true":
        return True
    kinds = {f.kind for f in rows}
    if {"start_codon", "stop_codon"} <= kinds:
        return True
    cds = [f for f in rows if f.kind == "CDS"]
    return bool(cds) and all(f.phase in ("0", "1", "2") for f in cds) and \
        sum(f.end - f.start + 1 for f in cds) % 3 == 0


def gene_models(gene):
    result = []
    tx_rows = {f.attr.get("ID"): f for f in gene.rows
               if f.kind in ("mRNA", "transcript")}
    for tid, children in gene.transcripts.items():
        rows = [tx_rows[tid]] + children
        introns = cds_introns(rows)
        if introns:
            result.append((tid, introns, rows))
    return result


def load_junctions(paths, min_unique, min_samples):
    samples, reads = defaultdict(set), Counter()
    for sample_no, path in enumerate(paths):
        with open(path) as handle:
            for lineno, line in enumerate(handle, 1):
                fields = line.split()
                if len(fields) < 9:
                    raise ValueError(f"{path}:{lineno}: malformed STAR junction row")
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                strand = {1: "+", 2: "-"}.get(int(fields[3]))
                motif, unique = int(fields[4]), int(fields[6])
                if strand and motif in range(1, 7) and unique >= min_unique:
                    key = (chrom, strand, start, end)
                    samples[key].add(sample_no)
                    reads[key] += unique
    return {key for key, seen in samples.items() if len(seen) >= min_samples}


def support_rate(gene, introns, junctions):
    return sum((gene.seqid, gene.strand, start, end) in junctions
               for start, end in introns) / len(introns)


def overlap(a, b):
    return a.seqid == b.seqid and a.strand == b.strand and \
        a.start <= b.end and b.start <= a.end


def rescue(evm_genes, track_genes, junctions, min_rate_delta=.25,
           max_exon_factor=1.5, max_exon_add=5):
    reasons, replacements = Counter(), {}
    for candidate in track_genes:
        hits = [gene for gene in evm_genes if overlap(candidate, gene)]
        if len(hits) != 1:
            reasons["fusion_or_no_unique_locus"] += 1
            continue
        evm = hits[0]
        c_models, e_models = gene_models(candidate), gene_models(evm)
        if not c_models or not e_models:
            reasons["no_cds_intron_chain"] += 1
            continue
        complete = [m for m in c_models if complete_orf(m[2])]
        if not complete:
            reasons["incomplete_orf"] += 1
            continue
        c_best = max(complete, key=lambda m: (support_rate(candidate, m[1], junctions), -len(m[1]), m[0]))
        e_best = max(e_models, key=lambda m: (support_rate(evm, m[1], junctions), -len(m[1]), m[0]))
        if c_best[1] == e_best[1]:
            reasons["same_chain"] += 1
            continue
        c_exons, e_exons = len(c_best[1]) + 1, len(e_best[1]) + 1
        if c_exons > e_exons * max_exon_factor or c_exons > e_exons + max_exon_add:
            reasons["exon_growth"] += 1
            continue
        delta = support_rate(candidate, c_best[1], junctions) - support_rate(evm, e_best[1], junctions)
        if delta < min_rate_delta:
            reasons["insufficient_support_delta"] += 1
            continue
        previous = replacements.get(evm.gene_id)
        rank = (support_rate(candidate, c_best[1], junctions), -c_exons,
                candidate.gene_id)
        if previous is None or rank > previous[0]:
            if previous is not None:
                reasons["superseded_candidate"] += 1
            replacements[evm.gene_id] = (rank, candidate)
        else:
            reasons["superseded_candidate"] += 1
    return {gid: item[1] for gid, item in replacements.items()}, reasons


def write_result(headers, evm_genes, replacements, output):
    outdir = os.path.dirname(os.path.abspath(output)) or "."
    fd, temporary = tempfile.mkstemp(dir=outdir, suffix=".rescue.tmp")
    try:
        with os.fdopen(fd, "w") as handle:
            for line in headers:
                handle.write(line + "\n")
            for gene in evm_genes:
                chosen = replacements.get(gene.gene_id, gene)
                for feat in chosen.rows:
                    handle.write(feat.raw + "\n")
        os.replace(temporary, output)
    except BaseException:
        if os.path.exists(temporary):
            os.unlink(temporary)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("evm_gff3")
    parser.add_argument("output_gff3")
    parser.add_argument("--track", action="append", required=True,
                        help="single-track GFF3; repeat for Liftoff/Helixer/AUGUSTUS")
    parser.add_argument("--star-sj", action="append", required=True,
                        help="STAR SJ.out.tab; repeat for samples")
    parser.add_argument("--min-unique", type=int, default=3)
    parser.add_argument("--min-samples", type=int, default=2)
    parser.add_argument("--min-rate-delta", type=float, default=.25)
    parser.add_argument("--max-exon-factor", type=float, default=1.5)
    parser.add_argument("--max-exon-add", type=int, default=5)
    args = parser.parse_args(argv)
    if not 0 < args.min_rate_delta <= 1 or args.min_unique < 1 or args.min_samples < 1:
        parser.error("support thresholds must be positive and rate delta <= 1")

    headers, evm = load_gff(args.evm_gff3)
    tracks = []
    for path in args.track:
        tracks.extend(load_gff(path)[1])
    junctions = load_junctions(args.star_sj, args.min_unique, args.min_samples)
    replacements, reasons = rescue(evm, tracks, junctions,
                                   args.min_rate_delta,
                                   args.max_exon_factor,
                                   args.max_exon_add)
    write_result(headers, evm, replacements, args.output_gff3)
    print(f"rescue_chains: replaced={len(replacements)} not_replaced={sum(reasons.values())}", file=sys.stderr)
    for reason, count in sorted(reasons.items()):
        print(f"  not_replaced {reason}={count}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
