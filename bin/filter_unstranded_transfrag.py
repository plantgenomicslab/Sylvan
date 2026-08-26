#!/usr/bin/env python3
"""Filter unstranded (single-exon) transfrag ORF models before they reach EVM.

transfragDump splits its output by strand and everything downstream used to read
only the stranded half. The unstranded half is not junk -- it is the single-exon
transcripts, which transfrag cannot orient because it infers strand from splice
junctions. Arabidopsis has 6,669 strictly single-exon genes, so discarding this
half wholesale loses real signal.

It is, however, the noisier half: intronic fragments, intergenic transcription and
unspliced pre-mRNA all land here. This filter keeps models that look like genes
rather than fragments, using only evidence already present in the GFF3 -- no
reference is consulted.
"""
import argparse
import re
import sys
from collections import defaultdict


def get_attr(text, key):
    m = re.search(r"(?:^|;)" + re.escape(key) + r"=([^;]+)", text)
    return m.group(1) if m else None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("gff3")
    ap.add_argument("--min-aa", type=int, default=150,
                    help="minimum CDS length in codons (default: 150)")
    ap.add_argument("--min-fragment-count", type=int, default=50,
                    help="minimum transfrag FragmentCount (default: 50)")
    ap.add_argument("--require-complete", action="store_true", default=True,
                    help="keep only Integrity=complete models (default: on)")
    ap.add_argument("--allow-partial", dest="require_complete",
                    action="store_false",
                    help="also keep models without Integrity=complete")
    args = ap.parse_args()

    records = {}
    order = []
    cds_length = defaultdict(int)
    gene_of_mrna = {}
    current_gene = None

    with open(args.gff3) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            if fields[2] == "gene":
                current_gene = get_attr(fields[8], "ID")
                if not current_gene:
                    raise SystemExit(f"gene row without ID: {line.rstrip()}")
                order.append(current_gene)
                records[current_gene] = [line]
                continue

            parent = get_attr(fields[8], "Parent")
            if not parent:
                continue

            if fields[2] == "mRNA":
                if parent in records:
                    gene_of_mrna[get_attr(fields[8], "ID")] = parent
                    records[parent].append(line)
                continue

            gene_id = gene_of_mrna.get(parent, parent if parent in records else None)
            if gene_id is None:
                continue
            records[gene_id].append(line)
            if fields[2] == "CDS":
                cds_length[gene_id] += int(fields[4]) - int(fields[3]) + 1

    kept = 0
    dropped = defaultdict(int)
    for gene_id in order:
        gene_line = records[gene_id][0]
        attrs = gene_line.rstrip("\n").split("\t")[8]

        if args.require_complete and get_attr(attrs, "Integrity") != "complete":
            dropped["not_complete"] += 1
            continue
        if cds_length[gene_id] < args.min_aa * 3:
            dropped["short_cds"] += 1
            continue
        try:
            fragment_count = int(get_attr(attrs, "FragmentCount") or 0)
        except ValueError:
            fragment_count = 0
        if fragment_count < args.min_fragment_count:
            dropped["low_support"] += 1
            continue

        sys.stdout.writelines(records[gene_id])
        kept += 1

    summary = " ".join(f"{k}={v}" for k, v in sorted(dropped.items()))
    print(f"filter_unstranded_transfrag: input={len(order)} kept={kept} "
          f"min_aa={args.min_aa} min_fragment_count={args.min_fragment_count} "
          f"dropped[{summary}]", file=sys.stderr)

    if kept == 0:
        raise SystemExit("filter retained zero gene models -- refusing to "
                         "hand an empty set downstream")


if __name__ == "__main__":
    main()
