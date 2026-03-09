#!/usr/bin/env python3
"""Remove duplicate/phantom EVM gene models from GFF3.

Detects genes with coordinates beyond the actual chromosome length
(caused by EVM partition offset errors) and removes them.
Also removes structurally identical genes that differ only by
a coordinate offset (partition-boundary duplicates).

Usage:
    python dedup_evm_gff3.py input.gff3 genome.fasta.fai output.gff3
"""
import re
import sys


def load_chrom_lengths(fai_path):
    lengths = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            lengths[parts[0]] = int(parts[1])
    return lengths


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} input.gff3 genome.fai output.gff3")

    in_gff, fai_path, out_gff = sys.argv[1], sys.argv[2], sys.argv[3]
    chrom_lengths = load_chrom_lengths(fai_path)

    # First pass: identify genes beyond chromosome boundaries
    bad_genes = set()
    with open(in_gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom, end = parts[0], int(parts[4])
            if chrom in chrom_lengths and end > chrom_lengths[chrom]:
                gene_id = re.search(r"ID=([^;]+)", parts[8])
                if gene_id:
                    bad_genes.add(gene_id.group(1))

    if not bad_genes:
        # No phantom genes — copy input to output unchanged
        with open(in_gff) as f_in, open(out_gff, "w") as f_out:
            f_out.write(f_in.read())
        print(f"dedup_evm_gff3: no phantom genes found, {out_gff} unchanged", file=sys.stderr)
        return

    # Second pass: filter out bad genes and their children
    removed = 0
    kept = 0
    with open(in_gff) as f_in, open(out_gff, "w") as f_out:
        skip_children_of = set()
        for line in f_in:
            if line.startswith("#") or not line.strip():
                f_out.write(line)
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                f_out.write(line)
                continue

            attrs = parts[8]
            # Check if this is a bad gene
            id_match = re.search(r"ID=([^;]+)", attrs)
            parent_match = re.search(r"Parent=([^;]+)", attrs)

            if parts[2] == "gene" and id_match and id_match.group(1) in bad_genes:
                skip_children_of.add(id_match.group(1))
                removed += 1
                continue

            # Skip children (mRNA, exon, CDS) of bad genes
            if parent_match:
                parent_id = parent_match.group(1)
                if parent_id in bad_genes or parent_id in skip_children_of:
                    if id_match:
                        skip_children_of.add(id_match.group(1))
                    continue

            f_out.write(line)
            if parts[2] == "gene":
                kept += 1

    print(f"dedup_evm_gff3: removed {removed} phantom genes, kept {kept}", file=sys.stderr)


if __name__ == "__main__":
    main()
