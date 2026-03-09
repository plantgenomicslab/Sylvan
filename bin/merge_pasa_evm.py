#!/usr/bin/env python3
"""Merge PASA-updated gene models with EVM genes that PASA didn't touch.

PASA annotation comparison only outputs genes with transcript evidence overlap.
EVM genes without transcript support are silently dropped. This script merges
the PASA output with the original EVM genes, keeping PASA-updated versions
where available and adding back EVM genes that PASA missed.

When PASA has structurally altered an EVM gene (chimeric merge or significant
CDS change), both proteins are compared via DIAMOND against SwissProt (reviewed).
Best hit is determined by lowest E-value, with alignment length as tie-breaker.

- If PASA wins or ties: keep PASA gene (preserves isoforms/UTRs)
- If EVM wins: use EVM as the gene, graft PASA isoforms that fall within
  EVM gene boundaries as additional mRNA isoforms
- If no overlap with PASA: add EVM gene as-is (rescued)

Usage:
    merge_pasa_evm.py pasa.gff3 evm.gff3 output.gff3 [genome.fa swissprot.fa]
"""
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict


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
                gene_id = re.search(r"ID=([^;]+)", parts[8])
                gid = gene_id.group(1) if gene_id else f"unknown_{len(genes)}"
                current_gene = (parts[0], int(parts[3]), int(parts[4]), parts[6], gid)
                current_lines = [line]
            elif current_gene:
                current_lines.append(line)

    if current_gene:
        genes.append(current_gene + (current_lines,))

    return genes


def count_cds(lines):
    """Count CDS features in a gene's GFF3 lines."""
    count = 0
    for line in lines:
        parts = line.strip().split("\t")
        if len(parts) >= 3 and parts[2] == "CDS":
            count += 1
    return count


def get_mrna_ids(lines):
    """Get mRNA IDs from gene's GFF3 lines."""
    mrna_ids = []
    for line in lines:
        parts = line.strip().split("\t")
        if len(parts) >= 9 and parts[2] == "mRNA":
            m = re.search(r"ID=([^;]+)", parts[8])
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
            mrna_id = re.search(r"ID=([^;]+)", parts[8])
            mid = mrna_id.group(1) if mrna_id else "unknown"
            current_mrna = (int(parts[3]), int(parts[4]), mid)
            current_mrna_lines = [line]
        elif current_mrna:
            current_mrna_lines.append(line)

    if current_mrna:
        mrnas.append(current_mrna + (current_mrna_lines,))

    return mrnas


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
        subprocess.run(
            ["gffread", gff_tmp, "-g", genome_fa, "-y", pep_tmp],
            capture_output=True, timeout=120
        )
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return {}

    if not os.path.exists(pep_tmp) or os.path.getsize(pep_tmp) == 0:
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


def build_evm_gene_with_pasa_isoforms(evm_lines, pasa_lines, evm_start, evm_end, evm_gene_id):
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
    for mrna_start, mrna_end, mrna_id, mrna_lines in pasa_mrnas:
        # Check if this PASA isoform is contained within EVM gene range
        if mrna_start >= evm_start and mrna_end <= evm_end:
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

    return output_lines, grafted


def main():
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} pasa.gff3 evm.gff3 output.gff3 [genome.fa swissprot.fa]",
              file=sys.stderr)
        sys.exit(1)

    pasa_gff, evm_gff, out_gff = sys.argv[1], sys.argv[2], sys.argv[3]
    genome_fa = sys.argv[4] if len(sys.argv) > 4 else None
    swissprot_fa = sys.argv[5] if len(sys.argv) > 5 else None

    use_diamond = genome_fa and swissprot_fa

    pasa_genes = parse_genes(pasa_gff)
    evm_genes = parse_genes(evm_gff)

    print(f"Loaded {len(pasa_genes)} PASA genes, {len(evm_genes)} EVM genes", file=sys.stderr)

    # Index PASA genes by (chrom, strand)
    pasa_by_loc = defaultdict(list)
    for i, (chrom, start, end, strand, gid, lines) in enumerate(pasa_genes):
        pasa_by_loc[(chrom, strand)].append((start, end, gid, lines, i))

    # First pass: classify each EVM gene
    rescued_lines = []  # EVM genes with no PASA coverage
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
            rescued_lines.extend(lines)
            continue

        # Check for structural alteration
        ps, pe, pgid, plines, pidx = best_pasa
        pasa_len = pe - ps + 1
        is_chimeric = pasa_len > 2 * evm_len
        evm_cds = count_cds(lines)
        pasa_cds = count_cds(plines)
        cds_changed = evm_cds > 0 and abs(pasa_cds - evm_cds) > max(1, evm_cds * 0.3)

        if is_chimeric or cds_changed:
            conflicts.append((ei, (chrom, start, end, strand, gid, lines),
                              pidx, (ps, pe, pgid, plines)))

    print(f"Found {rescued_count} EVM genes without PASA coverage (rescued)",
          file=sys.stderr)
    print(f"Found {len(conflicts)} structural conflicts to resolve", file=sys.stderr)

    # Second pass: resolve conflicts via DIAMOND
    # pasa_replaced[pidx] = merged_lines  — PASA gene replaced by EVM+PASA_isoforms
    pasa_replaced = {}
    stats = {"evm_wins": 0, "pasa_wins": 0, "tie": 0, "grafted_isoforms": 0}

    if conflicts and use_diamond:
        with tempfile.TemporaryDirectory() as tmpdir:
            diamond_db = make_diamond_db(swissprot_fa, tmpdir)
            if not diamond_db:
                print("WARNING: Could not create DIAMOND DB", file=sys.stderr)
                use_diamond = False

            if use_diamond:
                print(f"DIAMOND DB built, comparing {len(conflicts)} conflicts...",
                      file=sys.stderr)

                # Batch extract proteins
                evm_gene_map = {}
                pasa_gene_map = {}
                for ci, (ei, evm_info, pidx, pasa_info) in enumerate(conflicts):
                    evm_gene_map[f"evm_{ci}"] = evm_info[5]
                    pasa_gene_map[f"pasa_{ci}"] = pasa_info[3]

                evm_proteins = batch_extract_proteins(evm_gene_map, genome_fa, tmpdir, "evm")
                pasa_proteins = batch_extract_proteins(pasa_gene_map, genome_fa, tmpdir, "pasa")

                # Write query FASTAs
                evm_query = os.path.join(tmpdir, "evm_queries.fa")
                pasa_query = os.path.join(tmpdir, "pasa_queries.fa")

                with open(evm_query, "w") as f:
                    for ci in range(len(conflicts)):
                        evm_mrnas = get_mrna_ids(conflicts[ci][1][5])
                        for mid in evm_mrnas:
                            if mid in evm_proteins:
                                f.write(f">evm_{ci}\n{evm_proteins[mid]}\n")
                                break

                with open(pasa_query, "w") as f:
                    for ci in range(len(conflicts)):
                        pasa_mrnas = get_mrna_ids(conflicts[ci][3][3])
                        for mid in pasa_mrnas:
                            if mid in pasa_proteins:
                                f.write(f">pasa_{ci}\n{pasa_proteins[mid]}\n")
                                break

                # Run DIAMOND
                evm_hits = diamond_search(evm_query, diamond_db, tmpdir, "evm")
                pasa_hits = diamond_search(pasa_query, diamond_db, tmpdir, "pasa")

                # Resolve each conflict
                for ci, (ei, evm_info, pidx, pasa_info) in enumerate(conflicts):
                    evm_hit = evm_hits.get(f"evm_{ci}")
                    pasa_hit = pasa_hits.get(f"pasa_{ci}")
                    winner = compare_diamond(evm_hit, pasa_hit)

                    if winner == "pasa":
                        # PASA clearly better — keep PASA
                        stats["pasa_wins"] += 1
                    else:
                        # EVM wins or tie — use EVM gene (preserves CDS accuracy),
                        # graft PASA isoforms within EVM range
                        if winner == "evm":
                            stats["evm_wins"] += 1
                        else:
                            stats["tie"] += 1
                        chrom, start, end, strand, gid, lines = evm_info
                        merged_lines, n_grafted = build_evm_gene_with_pasa_isoforms(
                            lines, pasa_info[3], start, end, gid)
                        pasa_replaced[pidx] = merged_lines
                        stats["grafted_isoforms"] += n_grafted

    elif conflicts and not use_diamond:
        print("No DIAMOND DB — keeping PASA for all conflicts", file=sys.stderr)
        stats["tie"] = len(conflicts)

    # Write merged output
    with open(out_gff, "w") as f:
        f.write("##gff-version 3\n")
        for i, (_, _, _, _, _, lines) in enumerate(pasa_genes):
            if i in pasa_replaced:
                # Write the EVM gene + grafted PASA isoforms
                for line in pasa_replaced[i]:
                    f.write(line)
            else:
                for line in lines:
                    f.write(line)
        if rescued_lines:
            f.write("# Rescued EVM genes without PASA coverage\n")
            for line in rescued_lines:
                f.write(line)

    total_genes = len(pasa_genes) + rescued_count
    print(f"\nResults:", file=sys.stderr)
    print(f"  PASA genes: {len(pasa_genes)} ({len(pasa_replaced)} replaced by EVM)",
          file=sys.stderr)
    print(f"  Rescued EVM genes: {rescued_count}", file=sys.stderr)
    print(f"  Total genes: {total_genes}", file=sys.stderr)
    print(f"  Conflicts: {len(conflicts)} — EVM wins: {stats['evm_wins']}, "
          f"PASA wins: {stats['pasa_wins']}, tie: {stats['tie']}", file=sys.stderr)
    print(f"  PASA isoforms grafted into EVM genes: {stats['grafted_isoforms']}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
