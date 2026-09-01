#!/usr/bin/env python3
"""Unit tests for junction-chain arbitration in merge_pasa_evm.py.

Covers the arbitration primitives, the graft gate, flag parsing, and one
end-to-end merge where the junction verdict must restore the EVM chain
without any DIAMOND database present (the pure-junction path).

Run: micromamba run -n sylvan python bin/test_merge_pasa_evm_junction.py
"""
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from merge_pasa_evm import (arbitrate_by_junctions, isoform_introns_supported,
                            junction_support, parse_args)

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  PASS  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


def _junctions(pairs, chrom="chr1", strand="+"):
    """Build the lookup structure load_splice_junctions() produces."""
    table = {}
    table[(chrom, strand)] = set(pairs)
    table[(chrom, ".")] = set(pairs)
    return table


def test_junction_support():
    junc = _junctions([(5501, 5999), (6501, 6999)])
    # Chain with both introns supported
    both = ((5000, 5500), (6000, 6500), (7000, 7500))
    check("support: fully supported chain", junction_support(both, "chr1", "+", junc) == (2, 2))
    # Chain with one novel intron
    half = ((5000, 5500), (6200, 6500), (7000, 7500))
    check("support: half supported chain", junction_support(half, "chr1", "+", junc) == (1, 2))
    # Single-exon chain carries no signal
    check("support: single exon is None", junction_support(((5000, 7500),), "chr1", "+", junc) is None)


def test_arbitrate():
    junc = _junctions([(5501, 5999), (6501, 6999)])
    evm = ((5000, 5500), (6000, 6500), (7000, 7500))       # orphans: 0
    pasa = ((5000, 5500), (6200, 6500), (7000, 7500))      # orphans: 1
    # Orphan-count criterion decides regardless of delta
    check("arbitrate: fewer orphans wins (EVM)", arbitrate_by_junctions(evm, pasa, "chr1", "+", junc, 0.25) == "evm")
    check("arbitrate: fewer orphans wins (PASA)", arbitrate_by_junctions(pasa, evm, "chr1", "+", junc, 0.25) == "pasa")
    check("arbitrate: count decides even at large delta",
          arbitrate_by_junctions(evm, pasa, "chr1", "+", junc, 0.6) == "evm")
    junc2 = _junctions([(5501, 5999)])
    evm_short = ((5000, 5500), (6000, 6500), (7000, 7500))   # 1/2 supported, 1 orphan
    pasa_long = ((5000, 5500), (6000, 6500), (6800, 6900), (7000, 7500))  # 1/3 supported, 2 orphans
    check("arbitrate: more orphans loses even if longer",
          arbitrate_by_junctions(evm_short, pasa_long, "chr1", "+", junc2, 0.25) == "evm")
    # Equal orphans (1 vs 1), fractions 1/2 vs 2/3: diff 0.17 < 0.25 -> None
    junc3 = _junctions([(5501, 5999), (6501, 6799)])
    pasa_eq = ((5000, 5500), (6000, 6500), (6800, 6900), (7000, 7500))  # introns (6501,6799),(6901,6999): 1 sup 1 orphan
    check("arbitrate: equal orphans within delta is None",
          arbitrate_by_junctions(evm_short, pasa_eq, "chr1", "+", junc3, 0.25) is None)
    # Single-exon side -> indeterminate
    single = ((5000, 7500),)
    check("arbitrate: single-exon side is None", arbitrate_by_junctions(evm, single, "chr1", "+", junc, 0.25) is None)


def test_graft_gate():
    junc = _junctions([(5501, 5999)])
    supported = [
        "chr1\tpasa\tmRNA\t5000\t6500\t.\t+\t.\tID=iso1;Parent=g\n",
        "chr1\tpasa\tCDS\t5000\t5500\t.\t+\t0\tID=iso1.c1;Parent=iso1\n",
        "chr1\tpasa\tCDS\t6000\t6500\t.\t+\t0\tID=iso1.c2;Parent=iso1\n",
    ]
    novel = [
        "chr1\tpasa\tmRNA\t5000\t6500\t.\t+\t.\tID=iso2;Parent=g\n",
        "chr1\tpasa\tCDS\t5000\t5400\t.\t+\t0\tID=iso2.c1;Parent=iso2\n",
        "chr1\tpasa\tCDS\t6100\t6500\t.\t+\t0\tID=iso2.c2;Parent=iso2\n",
    ]
    single = ["chr1\tpasa\tmRNA\t5000\t5500\t.\t+\t.\tID=iso3;Parent=g\n",
              "chr1\tpasa\tCDS\t5000\t5500\t.\t+\t0\tID=iso3.c1;Parent=iso3\n"]
    check("graft: supported isoform passes", isoform_introns_supported(supported, "chr1", "+", junc))
    check("graft: novel-intron isoform rejected", not isoform_introns_supported(novel, "chr1", "+", junc))
    check("graft: single-exon isoform passes", isoform_introns_supported(single, "chr1", "+", junc))


def test_parse_args():
    pos, splice, delta = parse_args(["a.gff3", "b.gff3", "c.gff3",
                                     "--splice", "sp.tsv", "g.fa", "sw.fa",
                                     "--junction-delta", "0.4"])
    check("args: positionals preserved in order", pos == ["a.gff3", "b.gff3", "c.gff3", "g.fa", "sw.fa"])
    check("args: splice captured", splice == "sp.tsv")
    check("args: delta captured", delta == 0.4)
    pos2, splice2, delta2 = parse_args(["a", "b", "c"])
    check("args: defaults", splice2 is None and delta2 == 0.25 and pos2 == ["a", "b", "c"])


def _gene_block(source, gid, chain, chrom="chr1", strand="+"):
    start, end = chain[0][0], chain[-1][1]
    lines = [f"{chrom}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gid}\n",
             f"{chrom}\t{source}\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={gid}.m1;Parent={gid}\n"]
    for i, (s, e) in enumerate(chain, 1):
        lines.append(f"{chrom}\t{source}\tCDS\t{s}\t{e}\t.\t{strand}\t0\tID={gid}.c{i};Parent={gid}.m1\n")
    return "".join(lines)


def test_end_to_end_junction_restores_evm():
    """Junction-decisive conflict must restore the EVM chain with no DIAMOND DB."""
    evm_chain = ((5000, 5500), (6000, 6500), (7000, 7500))
    pasa_chain = ((5000, 5500), (6200, 6500), (7000, 7500))
    with tempfile.TemporaryDirectory() as tmp:
        pasa_gff = os.path.join(tmp, "pasa.gff3")
        evm_gff = os.path.join(tmp, "evm.gff3")
        splice = os.path.join(tmp, "trusted_splice")
        out = os.path.join(tmp, "merged.gff3")
        with open(pasa_gff, "w", encoding="utf-8") as f:
            f.write(_gene_block("pasa", "gA", pasa_chain))
        with open(evm_gff, "w", encoding="utf-8") as f:
            f.write(_gene_block("evm", "gA", evm_chain))
        with open(splice, "w", encoding="utf-8") as f:
            # trusted_splice records flanking exon boundaries (exon_end,
            # next_exon_start); the loader converts to intron spans
            # (5501-5999, 6501-6999) at load time.
            f.write("chr1 5500 6000 10 + 10 0 0 0\n")
            f.write("chr1 6500 7000 10 + 10 0 0 0\n")
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "merge_pasa_evm.py")
        result = subprocess.run(
            [sys.executable, script, pasa_gff, evm_gff, out, "--splice", splice],
            capture_output=True, text=True, timeout=60, check=False)
        with open(out, encoding="utf-8") as fh:
            merged = fh.read()
        check("e2e: exits cleanly", result.returncode == 0, result.stderr[-300:])
        check("e2e: junction verdict logged", "Junction arbitration: 1 EVM" in result.stderr,
              result.stderr[-300:])
        check("e2e: EVM chain present (6000 CDS)", "\t6000\t6500\t" in merged)
        check("e2e: PASA novel intron gone (6200 CDS)", "\t6200\t6500\t" not in merged)


def test_loader_conversion_and_strand():
    """Loader converts flank-exon boundaries to intron spans, strand-pure."""
    from merge_pasa_evm import load_splice_junctions
    with tempfile.TemporaryDirectory() as tmp:
        sp = os.path.join(tmp, "ts")
        with open(sp, "w", encoding="utf-8") as f:
            f.write("chr1 5500 6000 10 + 10 0 0 0\n")
            f.write("chr1 8800 9200 10 . 10 0 0 0\n")
        j = load_splice_junctions(sp)
        check("loader: flank->intron conversion", (5501, 5999) in j[("chr1", "+")])
        check("loader: stranded row NOT in '.' bucket", (5501, 5999) not in j.get(("chr1", "."), set()))
        check("loader: unstranded row in '.' bucket", (8801, 9199) in j[("chr1", ".")])
        check("loader: unstranded row not in '+' bucket", (8801, 9199) not in j.get(("chr1", "+"), set()))


def _script():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "merge_pasa_evm.py")


def test_empty_splice_falls_back():
    """A splice file with no usable rows must disable arbitration (keep PASA)."""
    evm_chain = ((5000, 5500), (6000, 6500))
    pasa_chain = ((5000, 5500), (6200, 6500))
    with tempfile.TemporaryDirectory() as tmp:
        pasa_gff, evm_gff = os.path.join(tmp, "p.gff3"), os.path.join(tmp, "e.gff3")
        splice, out = os.path.join(tmp, "empty_splice"), os.path.join(tmp, "m.gff3")
        with open(pasa_gff, "w", encoding="utf-8") as f:
            f.write(_gene_block("pasa", "gA", pasa_chain))
        with open(evm_gff, "w", encoding="utf-8") as f:
            f.write(_gene_block("evm", "gA", evm_chain))
        open(splice, "w", encoding="utf-8").close()
        result = subprocess.run([sys.executable, _script(), pasa_gff, evm_gff, out,
                                 "--splice", splice],
                                capture_output=True, text=True, timeout=60, check=False)
        with open(out, encoding="utf-8") as fh:
            merged = fh.read()
        check("empty splice: warns and disables", "no usable junctions" in result.stderr,
              result.stderr[-300:])
        check("empty splice: PASA kept (legacy tie)", "\t6200\t6500\t" in merged)


def test_missing_flag_value_errors():
    result = subprocess.run([sys.executable, _script(), "a", "b", "c", "--splice"],
                            capture_output=True, text=True, timeout=30, check=False)
    check("args: missing --splice value exits 2", result.returncode == 2, str(result.returncode))


def test_mixed_verdict_chimera_backfill():
    """One PASA chimera vs two EVM genes with mixed verdicts keeps BOTH loci.

    gA: junction says EVM (PASA primary chain carries 2 orphans vs 0).
    gB: orphan counts tie, no DIAMOND -> verdict stays PASA/tie — but gA's win
    removes the chimera block, so gB's EVM model must be backfilled.
    """
    with tempfile.TemporaryDirectory() as tmp:
        pasa_gff, evm_gff = os.path.join(tmp, "p.gff3"), os.path.join(tmp, "e.gff3")
        splice, out = os.path.join(tmp, "ts"), os.path.join(tmp, "m.gff3")
        # PASA chimera: one gene 5000-19500, primary chain spans both loci
        chim = ["chr1\tpasa\tgene\t5000\t19500\t.\t+\t.\tID=chimera\n",
                "chr1\tpasa\tmRNA\t5000\t19500\t.\t+\t.\tID=chimera.m1;Parent=chimera\n"]
        for i, (s, e) in enumerate(((5000, 5500), (6200, 6500), (17000, 17400)), 1):
            chim.append(f"chr1\tpasa\tCDS\t{s}\t{e}\t.\t+\t0\tID=chimera.c{i};Parent=chimera.m1\n")
        with open(pasa_gff, "w", encoding="utf-8") as f:
            f.write("".join(chim))
        with open(evm_gff, "w", encoding="utf-8") as f:
            f.write(_gene_block("evm", "gA", ((5000, 5500), (6000, 6500))))
            f.write(_gene_block("evm", "gB", ((17000, 17400), (17800, 18200), (18800, 19200))))
        with open(splice, "w", encoding="utf-8") as f:
            f.write("chr1 5500 6000 10 + 10 0 0 0\n")   # supports gA's intron only
        result = subprocess.run([sys.executable, _script(), pasa_gff, evm_gff, out,
                                 "--splice", splice],
                                capture_output=True, text=True, timeout=60, check=False)
        with open(out, encoding="utf-8") as fh:
            merged = fh.read()
        check("chimera: exits cleanly", result.returncode == 0, result.stderr[-300:])
        check("chimera: gA EVM chain present", "\t6000\t6500\t" in merged)
        check("chimera: gB locus preserved (backfilled)", "\t17800\t18200\t" in merged)
        check("chimera: PASA novel intron gone", "\t6200\t6500\t" not in merged)
        check("chimera: backfill logged", "backfilled" in result.stderr, result.stderr[-400:])


def main():
    test_junction_support()
    test_arbitrate()
    test_graft_gate()
    test_parse_args()
    test_loader_conversion_and_strand()
    test_end_to_end_junction_restores_evm()
    test_empty_splice_falls_back()
    test_missing_flag_value_errors()
    test_mixed_verdict_chimera_backfill()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILED: {FAILURES}")
        sys.exit(1)
    print("\nAll tests passed")


if __name__ == "__main__":
    main()
