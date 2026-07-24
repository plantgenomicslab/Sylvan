#!/usr/bin/env python3
"""Drop coordinate-invalid gene models from an AUGUSTUS training GFF3.

Background (tracker wyim-pgl/Sylvan-EGAPx#31): the GETA training set fed to
BGM2AT can contain a handful of models whose CDS exons butt against or
overlap each other, giving introns of length <= 0. AUGUSTUS's GBProcessor
skips such models but logs one error per re-read — with 5 rounds x 12-fold
cross-validation the same 3-4 genes produce ~9k log lines, masking real
errors.

An mRNA is invalid when, after sorting its CDS features by genomic start,
any adjacent pair has ``next.start <= prev.end`` (intron <= 0 bp). Invalid
mRNAs are dropped together with their child features; a gene row is dropped
only when no mRNA remains under it.

Usage:
    python sanitize_training_gff.py input.gff3 output.gff3
"""
import re
import sys

CDS_PARENT_RE = re.compile(r"Parent=([^;\s]+)")
ID_RE = re.compile(r"ID=([^;\s]+)")


def _attr(regex, attrs):
    match = regex.search(attrs)
    return match.group(1) if match else None


def find_invalid_mrnas(lines):
    """Return the set of mRNA IDs whose CDS layout has an intron <= 0 bp."""
    mrna_cds = {}
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "CDS":
            continue
        parent = _attr(CDS_PARENT_RE, fields[8])
        if parent is None:
            continue
        start, end = int(fields[3]), int(fields[4])
        mrna_cds.setdefault(parent, []).append((start, end))

    invalid = set()
    for mrna, cds_list in mrna_cds.items():
        cds_list.sort()
        for (prev_start, prev_end), (next_start, _next_end) in zip(cds_list, cds_list[1:]):
            if next_start <= prev_end:
                invalid.add(mrna)
                break
    return invalid


def sanitize_gff(lines):
    """Remove invalid mRNAs (and orphaned gene rows) from GFF3 lines.

    Returns (cleaned_lines, dropped_mrna_ids).
    """
    invalid = find_invalid_mrnas(lines)
    if not invalid:
        return list(lines), set()

    gene_to_mrnas = {}
    kept = []
    for line in lines:
        if line.startswith("#"):
            kept.append(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            kept.append(line)
            continue
        ftype, attrs = fields[2], fields[8]
        row_id = _attr(ID_RE, attrs)
        parent = _attr(CDS_PARENT_RE, attrs)

        if ftype == "gene":
            kept.append(line)  # filtered below once gene_to_mrnas is complete
            continue
        if ftype == "mRNA":
            if row_id in invalid:
                continue
            if parent:
                gene_to_mrnas.setdefault(parent, set()).add(row_id)
            kept.append(line)
            continue
        # Child features (CDS/exon/UTR/...) — drop if their parent is invalid.
        if parent in invalid:
            continue
        kept.append(line)

    # Gene rows: keep only those that still have >= 1 mRNA in the output.
    # mRNAs may carry their Parent= even when valid, so map ALL mRNA rows
    # (not just kept ones) then subtract invalid.
    all_gene_mrnas = {}
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "mRNA":
            continue
        parent = _attr(CDS_PARENT_RE, fields[8])
        row_id = _attr(ID_RE, fields[8])
        if parent and row_id:
            all_gene_mrnas.setdefault(parent, set()).add(row_id)

    final = []
    for line in kept:
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 9 and fields[2] == "gene":
            gene_id = _attr(ID_RE, fields[8])
            mrnas = all_gene_mrnas.get(gene_id, set())
            if mrnas and mrnas <= invalid:
                continue  # every isoform was dropped
        final.append(line)
    return final, invalid


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} input.gff3 output.gff3")
    with open(sys.argv[1]) as fh:
        lines = fh.read().splitlines()
    cleaned, dropped = sanitize_gff(lines)
    with open(sys.argv[2], "w") as fh:
        fh.write("\n".join(cleaned) + "\n")
    print(
        f"sanitize_training_gff: kept {len(lines) - len(cleaned)} invalid-feature lines; "
        f"dropped {len(dropped)} mRNAs: {sorted(dropped)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
