"""Unit tests for merge_evm_partitions.py.

Synthetic layout: two segmented partitions of T_Chr1 with the production
20 kb overlap (1-1000000 and 980001-1980000) plus one unsegmented
chromosome (T_ChrC, 3-field partitions_list line). Coordinates inside a
partition file are partition-local, exactly what EVM_to_GFF3.pl emits.
"""
import os
import subprocess
import sys
import tempfile
import unittest

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "merge_evm_partitions.py")


def gff_model(chrom, gid, mid, spans, strand="+"):
    """One EVM-style gene/mRNA with paired exon+CDS rows per span."""
    lo = min(s for s, _ in spans)
    hi = max(e for _, e in spans)
    rows = [
        f"{chrom}\tEVM\tgene\t{lo}\t{hi}\t.\t{strand}\t.\tID={gid}",
        f"{chrom}\tEVM\tmRNA\t{lo}\t{hi}\t.\t{strand}\t.\tID={mid};Parent={gid}",
    ]
    for i, (s, e) in enumerate(spans, 1):
        rows.append(f"{chrom}\tEVM\texon\t{s}\t{e}\t.\t{strand}\t."
                    f"\tID={mid}.exon{i};Parent={mid}")
        rows.append(f"{chrom}\tEVM\tCDS\t{s}\t{e}\t.\t{strand}\t0"
                    f"\tID=cds.{mid};Parent={mid}")
    return "\n".join(rows) + "\n"


class TestMergeEvmPartitions(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = self._tmp.name
        self.p1 = os.path.join(self.root, "T_Chr1", "T_Chr1_1-1000000")
        self.p2 = os.path.join(self.root, "T_Chr1", "T_Chr1_980001-1980000")
        self.pc = os.path.join(self.root, "T_ChrC")
        for d in (self.p1, self.p2, self.pc):
            os.makedirs(d)
        with open(os.path.join(self.root, "partitions_list.out"), "w") as fh:
            fh.write(f"T_Chr1\t{os.path.dirname(self.p1)}\tY\t{self.p1}\n")
            fh.write(f"T_Chr1\t{os.path.dirname(self.p2)}\tY\t{self.p2}\n")
            fh.write(f"T_ChrC\t{self.pc}\tN\n")
        # every partition starts empty; tests overwrite what they need
        for d in (self.p1, self.p2, self.pc):
            self.write(d, "")

    def tearDown(self):
        self._tmp.cleanup()

    def write(self, part, text):
        with open(os.path.join(part, "evm.partition.gff3"), "w") as fh:
            fh.write(text)

    def run_merge(self):
        out = os.path.join(self.root, "merged.gff3")
        subprocess.run(
            [sys.executable, SCRIPT,
             os.path.join(self.root, "partitions_list.out"), out],
            check=True, capture_output=True, text=True)
        with open(out) as fh:
            return fh.read()

    def test_core_ownership_and_offset(self):
        # p1 core = [1, 980000]: local 5000-6000, mid 5500 -> p1 owns it.
        self.write(self.p1, gff_model("T_Chr1", "g1", "m1", [(5000, 6000)]))
        # p2 local 21000-22000 -> genome 1001000-1002000, inside p2 core.
        self.write(self.p2, gff_model("T_Chr1", "g2", "m2", [(21000, 22000)]))
        out = self.run_merge()
        self.assertIn("ID=T_Chr1_1-1000000_g1", out)
        self.assertIn("\t1001000\t1002000\t", out)  # offset 980000 applied
        self.assertIn("ID=T_Chr1_980001-1980000_g2", out)
        self.assertEqual(out.count("\tmRNA\t"), 2)

    def test_unsegmented_chromosome_kept_whole(self):
        self.write(self.pc, gff_model("T_ChrC", "gc", "mc", [(100, 400)]))
        out = self.run_merge()
        self.assertIn("ID=T_ChrC_gc", out)
        self.assertIn("\t100\t400\t", out)  # no offset

    def test_boundary_model_kept_once(self):
        # Genome 989000-991000, mid 990000: p2's core. Both partitions carry
        # a copy (each in its own local coordinates) -- exactly the state
        # stock recombine tears apart row by row.
        self.write(self.p1, gff_model("T_Chr1", "g1", "m1", [(989000, 991000)]))
        self.write(self.p2, gff_model("T_Chr1", "g1", "m1", [(9000, 11000)]))
        out = self.run_merge()
        self.assertEqual(out.count("\tmRNA\t"), 1)
        self.assertIn("T_Chr1_980001-1980000_m1", out)

    def test_v2_prefers_more_complete_neighbor(self):
        # Owner is p2: its 2-span copy sits at genome 990000-996000
        # (p2 local 10000-11000, 15000-16000; mid 993000 in p2 core).
        self.write(self.p2, gff_model("T_Chr1", "g1", "m1",
                                      [(10000, 11000), (15000, 16000)]))
        # p1 holds a more complete 4-span version of the same locus
        # (genome==local for p1). Its own mid 993000 is outside p1's core,
        # so it is not independently kept -- only the v2 swap can save it.
        self.write(self.p1, gff_model("T_Chr1", "g1", "m1full",
                                      [(989500, 990500), (992000, 992500),
                                       (994000, 994500), (995500, 996500)]))
        out = self.run_merge()
        self.assertEqual(out.count("\tmRNA\t"), 1)
        self.assertIn("T_Chr1_1-1000000_m1full", out)
        self.assertEqual(out.count("\tCDS\t"), 4)
        self.assertNotIn("_m1;", out)  # 2-span copy replaced

    def test_gene_row_emitted_once(self):
        self.write(self.p1, gff_model("T_Chr1", "g1", "m1", [(5000, 6000)]))
        out = self.run_merge()
        self.assertEqual(out.count("\tgene\t"), 1)

    def test_deterministic_under_row_shuffle(self):
        text = (gff_model("T_Chr1", "g1", "m1", [(5000, 6000)])
                + gff_model("T_Chr1", "g2", "m2",
                            [(8000, 9500), (9800, 9900)]))
        self.write(self.p1, text)
        out_a = self.run_merge()
        lines = [l for l in text.splitlines() if l]
        self.write(self.p1, "\n".join(reversed(lines)) + "\n")
        out_b = self.run_merge()
        self.assertEqual(out_a, out_b)
        self.assertEqual(out_a.count("\tmRNA\t"), 2)

    def test_missing_partition_file_is_tolerated(self):
        os.unlink(os.path.join(self.pc, "evm.partition.gff3"))
        self.write(self.p1, gff_model("T_Chr1", "g1", "m1", [(5000, 6000)]))
        out = self.run_merge()
        self.assertEqual(out.count("\tmRNA\t"), 1)


if __name__ == "__main__":
    unittest.main()
