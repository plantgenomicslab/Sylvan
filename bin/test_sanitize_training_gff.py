#!/usr/bin/env python3
"""Unit test for bin/sanitize_training_gff.py.

The GETA → AUGUSTUS training genbank (genes.gb) can contain a handful of
coordinate-invalid gene models whose CDS exons butt against or overlap each
other (intron length <= 0). AUGUSTUS's GBProcessor skips them but logs
~9k repeated errors per run (tracker issue wyim-pgl/Sylvan-EGAPx#31).

sanitize_training_gff.py must drop exactly those invalid mRNAs (and their
now-empty genes) while preserving every valid model byte-for-byte in
structure.
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sanitize_training_gff import (  # noqa: E402
    find_invalid_mrnas,
    sanitize_gff,
)


def _gff(chrom, source, ftype, start, end, strand, attrs):
    return "\t".join(
        [chrom, source, ftype, str(start), str(end), ".", strand, ".", attrs]
    )


def _model_lines(gene_id, tx_id, chrom, strand, cds_list, exon_list):
    """Build one gene+mRNA+exons+CDS model as GFF3 lines."""
    g_start = min(s for s, _ in exon_list)
    g_end = max(e for _, e in exon_list)
    lines = [
        _gff(chrom, ".", "gene", g_start, g_end, strand, f"ID={gene_id};"),
        _gff(
            chrom, ".", "mRNA", g_start, g_end, strand,
            f"ID={tx_id};Parent={gene_id};Integrity=complete;",
        ),
    ]
    for i, (s, e) in enumerate(exon_list, 1):
        lines.append(
            _gff(chrom, ".", "exon", s, e, strand,
                 f"ID={tx_id}.exon{i};Parent={tx_id};")
        )
    for i, (s, e) in enumerate(cds_list, 1):
        lines.append(
            _gff(chrom, ".", "CDS", s, e, strand,
                 f"ID={tx_id}.CDS{i};Parent={tx_id};")
        )
    return lines


class TestFindInvalidMrnas(unittest.TestCase):
    def test_zero_length_intron_minus_strand_real_case(self):
        """The real Ath genewise08931 case: CDS2 ends at 12427955 and CDS3
        also ends at 12427955 on the minus strand → intron length -1 → invalid."""
        lines = _model_lines(
            "genewise08931", "genewise08931.t1", "Ath_Chr2", "-",
            cds_list=[(12428249, 12428678), (12427955, 12427959), (12427689, 12427955)],
            exon_list=[(12428249, 12428678), (12427689, 12427959)],
        )
        invalid = find_invalid_mrnas(lines)
        self.assertEqual(invalid, {"genewise08931.t1"})

    def test_normal_gene_is_valid(self):
        """CDS separated by a real intron must NOT be flagged."""
        lines = _model_lines(
            "geneA", "geneA.t1", "Chr1", "+",
            cds_list=[(1000, 1200), (1500, 1700), (2000, 2200)],
            exon_list=[(1000, 1200), (1500, 1700), (2000, 2200)],
        )
        self.assertEqual(find_invalid_mrnas(lines), set())

    def test_touching_cds_plus_strand_invalid(self):
        """Adjacent CDS with 0-length intron on plus strand → invalid."""
        lines = _model_lines(
            "geneB", "geneB.t1", "Chr1", "+",
            cds_list=[(1000, 1246), (1246, 1250)],
            exon_list=[(1000, 1250)],
        )
        self.assertEqual(find_invalid_mrnas(lines), {"geneB.t1"})

    def test_overlapping_cds_invalid(self):
        """Overlapping CDS (negative intron) → invalid."""
        lines = _model_lines(
            "geneC", "geneC.t1", "Chr3", "+",
            cds_list=[(1000, 1500), (1400, 1800)],
            exon_list=[(1000, 1800)],
        )
        self.assertEqual(find_invalid_mrnas(lines), {"geneC.t1"})

    def test_one_bp_gap_is_valid(self):
        """A 1-bp intron is degenerate but structurally positive → kept."""
        lines = _model_lines(
            "geneD", "geneD.t1", "Chr1", "+",
            cds_list=[(1000, 1200), (1202, 1500)],
            exon_list=[(1000, 1500)],
        )
        self.assertEqual(find_invalid_mrnas(lines), set())

    def test_single_cds_gene_is_valid(self):
        lines = _model_lines(
            "geneE", "geneE.t1", "Chr5", "-",
            cds_list=[(500, 900)],
            exon_list=[(500, 900)],
        )
        self.assertEqual(find_invalid_mrnas(lines), set())


class TestSanitizeGff(unittest.TestCase):
    def setUp(self):
        self.bad = _model_lines(
            "genewise08931", "genewise08931.t1", "Ath_Chr2", "-",
            cds_list=[(12428249, 12428678), (12427955, 12427959), (12427689, 12427955)],
            exon_list=[(12428249, 12428678), (12427689, 12427959)],
        )
        self.good = _model_lines(
            "geneA", "geneA.t1", "Chr1", "+",
            cds_list=[(1000, 1200), (1500, 1700)],
            exon_list=[(1000, 1200), (1500, 1700)],
        )

    def test_invalid_gene_fully_removed_valid_preserved(self):
        """Invalid gene + mRNA + all child features gone; valid model intact."""
        lines = self.good + self.bad
        cleaned, dropped = sanitize_gff(lines)
        self.assertEqual(dropped, {"genewise08931.t1"})
        self.assertEqual(cleaned, self.good)

    def test_gene_row_dropped_only_when_no_mrna_remains(self):
        """A gene with one bad and one good isoform keeps the gene row."""
        good_iso = _model_lines(
            "geneF", "geneF.t2", "Chr2", "+",
            cds_list=[(5000, 5200), (5500, 5700)],
            exon_list=[(5000, 5200), (5500, 5700)],
        )
        # Rebuild geneF as a two-isoform gene: t1 (bad) + t2 (good)
        bad_iso = [
            _gff("Chr2", ".", "mRNA", 100, 300, "+",
                 "ID=geneF.t1;Parent=geneF;"),
            _gff("Chr2", ".", "CDS", 100, 200, "+",
                 "ID=geneF.t1.CDS1;Parent=geneF.t1;"),
            _gff("Chr2", ".", "CDS", 200, 300, "+",
                 "ID=geneF.t1.CDS2;Parent=geneF.t1;"),
        ]
        gene_row = _gff("Chr2", ".", "gene", 100, 5700, "+", "ID=geneF;")
        good_iso_no_gene = good_iso[1:]  # strip duplicate gene row
        lines = [gene_row] + bad_iso + good_iso_no_gene
        cleaned, dropped = sanitize_gff(lines)
        self.assertEqual(dropped, {"geneF.t1"})
        self.assertIn(gene_row, cleaned)
        for l in bad_iso:
            self.assertNotIn(l, cleaned)
        for l in good_iso_no_gene:
            self.assertIn(l, cleaned)

    def test_non_gff_lines_passthrough(self):
        """Comments/headers survive untouched."""
        lines = ["##gff-version 3", "# a comment"] + self.good
        cleaned, _ = sanitize_gff(lines)
        self.assertEqual(cleaned[0], "##gff-version 3")
        self.assertEqual(cleaned[1], "# a comment")


if __name__ == "__main__":
    unittest.main()
