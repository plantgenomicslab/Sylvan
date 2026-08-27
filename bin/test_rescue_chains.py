#!/usr/bin/env python3
"""Tests for exact-junction whole-locus rescue and its guardrails."""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rescue_chains import load_gff, load_junctions, rescue, write_result  # noqa: E402


def bundle(gid, starts, source="EVM", complete=True, seqid="Chr1"):
    end = starts[-1] + 29
    attr = "valid_ORF=True" if complete else "valid_ORF=False;missing_stop_codon=True"
    rows = [f"{seqid}\t{source}\tgene\t{starts[0]}\t{end}\t.\t+\t.\tID={gid}",
            f"{seqid}\t{source}\tmRNA\t{starts[0]}\t{end}\t.\t+\t.\tID={gid}.t1;Parent={gid};{attr}"]
    for n, start in enumerate(starts, 1):
        rows.append(f"{seqid}\t{source}\tCDS\t{start}\t{start+29}\t.\t+\t0\tID={gid}.c{n};Parent={gid}.t1")
    return "\n".join(rows) + "\n"


class TestRescueChains(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def put(self, name, text):
        path = os.path.join(self.tmp.name, name)
        with open(path, "w") as handle:
            handle.write(text)
        return path

    def genes(self, name, text):
        return load_gff(self.put(name, text))[1]

    def test_exact_supported_chain_replaces_whole_bundle(self):
        evm = self.genes("evm.gff", bundle("E", [100, 200]))
        cand = self.genes("cand.gff", bundle("H", [100, 250], "Helixer"))
        sj = self.put("SJ", "Chr1\t130\t249\t1\t1\t0\t5\t0\t20\n")
        junctions = load_junctions([sj], 3, 1)
        repl, reasons = rescue(evm, cand, junctions, min_rate_delta=.5)
        self.assertEqual(set(repl), {"E"})
        out = self.put("out.gff", "")
        write_result(["##gff-version 3"], evm, repl, out)
        with open(out) as handle:
            text = handle.read()
        self.assertIn("ID=H", text)
        self.assertNotIn("ID=E", text)
        self.assertFalse(reasons)

    def test_tolerance_zero_does_not_accept_near_junction(self):
        evm = self.genes("e", bundle("E", [100, 200]))
        cand = self.genes("c", bundle("H", [100, 250]))
        repl, reasons = rescue(evm, cand, {("Chr1", "+", 131, 249)})
        self.assertFalse(repl)
        self.assertEqual(reasons["insufficient_support_delta"], 1)

    def test_incomplete_orf_blocked(self):
        evm = self.genes("e", bundle("E", [100, 200]))
        cand = self.genes("c", bundle("H", [100, 250], complete=False))
        repl, reasons = rescue(evm, cand, {("Chr1", "+", 130, 249)})
        self.assertFalse(repl)
        self.assertEqual(reasons["incomplete_orf"], 1)

    def test_fusion_candidate_blocked(self):
        evm = self.genes("e", bundle("E1", [100, 200]) + bundle("E2", [260, 360]))
        cand = self.genes("c", bundle("H", [100, 360]))
        repl, reasons = rescue(evm, cand, {("Chr1", "+", 130, 359)})
        self.assertFalse(repl)
        self.assertEqual(reasons["fusion_or_no_unique_locus"], 1)

    def test_exon_explosion_blocked(self):
        evm = self.genes("e", bundle("E", [100, 500]))
        cand = self.genes("c", bundle("H", [100, 150, 200, 250, 300, 350, 400]))
        junctions = {("Chr1", "+", a + 30, b - 1)
                     for a, b in zip([100, 150, 200, 250, 300, 350],
                                     [150, 200, 250, 300, 350, 400])}
        repl, reasons = rescue(evm, cand, junctions)
        self.assertFalse(repl)
        self.assertEqual(reasons["exon_growth"], 1)


if __name__ == "__main__":
    unittest.main()
