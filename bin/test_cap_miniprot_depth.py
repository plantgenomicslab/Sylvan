#!/usr/bin/env python3
"""Unit test for bin/cap_miniprot_depth.py.

EVM 2.1.0 accumulates evidence weight per base (`add_match_coverage`) and per
exact intron key (`add_introns`) with no cap, average, or redundancy
normalisation. Feeding it the full OrthoDB Viridiplantae set therefore lets
miniprot dominate: a 1 Mb Ath_ODB partition carries 654,453 miniprot segments
from 124,450 distinct query proteins, at 297x mean depth over covered bases and
2,316x at the peak, against 1,075 GeneWise rows.

cap_miniprot_depth.py bounds both accumulation paths before EVM sees the file:

    per genomic base            retained coding segments <= k_base
    per exact intron key        supporting alignments    <= k_intron

Selection is atomic per alignment ID (one genomic hit, possibly many CDS
segments) -- never per Target (one query protein, possibly many hits) and never
per row, because Identity varies between segments of the same alignment in 99%
of multi-segment alignments.

Design: docs/superpowers/specs/2026-08-25-evm-selftune-design.md
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cap_miniprot_depth import (  # noqa: E402
    cap_miniprot_file,
    audit_selection,
    cap_miniprot_lines,
    select_alignments,
    summarize_alignments,
)


def _seg(aln_id, start, end, identity, seqid="C1", strand="+", rank=1,
         score=100, target="Q1", target_span=(1, 100)):
    """One miniprot CDS segment as an EVM-format GFF row."""
    t_start, t_end = target_span
    attrs = (f"ID={aln_id};Rank={rank};Identity={identity};"
             f"Target={target} {t_start} {t_end}")
    return "\t".join([seqid, "miniprot", "miniprot", str(start), str(end),
                      str(score), strand, "0", attrs])


def _spliced(aln_id, blocks, identity=0.9, **kw):
    """An alignment spanning several exon blocks -> introns between them."""
    return [_seg(aln_id, s, e, identity, **kw) for s, e in blocks]


def _depth_at(lines, seqid, strand, pos):
    """How many retained segments cover `pos`."""
    n = 0
    for ln in lines:
        f = ln.split("\t")
        if f[0] == seqid and f[6] == strand and int(f[3]) <= pos <= int(f[4]):
            n += 1
    return n


def _intron_support(lines):
    """Map exact intron key -> number of distinct alignments supporting it."""
    by_id = {}
    for ln in lines:
        f = ln.split("\t")
        aln = f[8].split("ID=")[1].split(";")[0]
        by_id.setdefault((f[0], f[6], aln), []).append((int(f[3]), int(f[4])))
    support = {}
    for (seqid, strand, aln), segs in by_id.items():
        segs.sort()
        for (_, e1), (s2, _) in zip(segs, segs[1:]):
            support.setdefault((seqid, strand, e1 + 1, s2 - 1), set()).add(aln)
    return {k: len(v) for k, v in support.items()}


def _ids(lines):
    return {ln.split("ID=")[1].split(";")[0] for ln in lines}


class TestPerBaseCap(unittest.TestCase):
    def test_depth_never_exceeds_k_base(self):
        """Three alignments stacked on one base, k_base=2 -> one is dropped."""
        lines = [_seg("MP1", 100, 200, 0.9),
                 _seg("MP2", 100, 200, 0.8),
                 _seg("MP3", 100, 200, 0.7)]

        kept = cap_miniprot_lines(lines, k_base=2, k_intron=99)

        self.assertEqual(_depth_at(kept, "C1", "+", 150), 2)

    def test_non_overlapping_alignments_all_survive(self):
        """A cap constrains overlap, not total count."""
        lines = [_seg("MP1", 100, 200, 0.9),
                 _seg("MP2", 300, 400, 0.9),
                 _seg("MP3", 500, 600, 0.9)]

        kept = cap_miniprot_lines(lines, k_base=1, k_intron=99)

        self.assertEqual(_ids(kept), {"MP1", "MP2", "MP3"})

    def test_opposite_strands_do_not_share_capacity(self):
        """EVM scores per (seqid, strand); capacity must follow."""
        lines = [_seg("MP1", 100, 200, 0.9, strand="+"),
                 _seg("MP2", 100, 200, 0.9, strand="-")]

        kept = cap_miniprot_lines(lines, k_base=1, k_intron=99)

        self.assertEqual(_ids(kept), {"MP1", "MP2"})


class TestIntronCap(unittest.TestCase):
    def test_intron_multiplicity_never_exceeds_k_intron(self):
        """Same exact intron supported 3x, k_intron=2 -> one alignment drops."""
        blocks = [(100, 200), (300, 400)]
        lines = (_spliced("MP1", blocks, identity=0.9)
                 + _spliced("MP2", blocks, identity=0.8)
                 + _spliced("MP3", blocks, identity=0.7))

        kept = cap_miniprot_lines(lines, k_base=99, k_intron=2)

        self.assertTrue(all(v <= 2 for v in _intron_support(kept).values()))

    def test_intron_cap_binds_when_base_cap_would_not(self):
        """The two constraints are not equivalent -- intron cap must bite alone."""
        blocks = [(100, 200), (300, 400)]
        lines = (_spliced("MP1", blocks, identity=0.9)
                 + _spliced("MP2", blocks, identity=0.8))

        kept = cap_miniprot_lines(lines, k_base=99, k_intron=1)

        self.assertEqual(len(_ids(kept)), 1)


class TestAtomicity(unittest.TestCase):
    def test_accepted_alignment_keeps_every_segment(self):
        """An admitted ID is admitted whole -- no chimeric partial alignment."""
        lines = _spliced("MP1", [(100, 200), (300, 400), (500, 600)])

        kept = cap_miniprot_lines(lines, k_base=5, k_intron=5)

        self.assertEqual(len(kept), 3)

    def test_rejected_alignment_leaves_no_segment(self):
        lines = (_spliced("MP1", [(100, 200), (300, 400)], identity=0.9)
                 + _spliced("MP2", [(100, 200), (300, 400)], identity=0.5))

        kept = cap_miniprot_lines(lines, k_base=1, k_intron=1)

        self.assertEqual(_ids(kept), {"MP1"})
        self.assertEqual(len(kept), 2)

    def test_same_target_different_alignments_are_independent(self):
        """One query protein hitting two loci must not be all-or-nothing."""
        lines = [_seg("MP1", 100, 200, 0.9, rank=1, target="Q1"),
                 _seg("MP2", 5000, 5100, 0.8, rank=2, target="Q1")]

        kept = cap_miniprot_lines(lines, k_base=1, k_intron=99)

        self.assertEqual(_ids(kept), {"MP1", "MP2"})


class TestQualityOrdering(unittest.TestCase):
    def test_higher_identity_wins_scarce_capacity(self):
        lines = [_seg("MP_low", 100, 200, 0.55),
                 _seg("MP_high", 100, 200, 0.95)]

        kept = cap_miniprot_lines(lines, k_base=1, k_intron=99)

        self.assertEqual(_ids(kept), {"MP_high"})

    def test_raw_score_does_not_outrank_identity(self):
        """miniprot score scales with alignment length; a long poor hit must
        not evict a short excellent one."""
        lines = [_seg("MP_long_poor", 100, 400, 0.55, score=9000),
                 _seg("MP_short_good", 100, 400, 0.95, score=100)]

        kept = cap_miniprot_lines(lines, k_base=1, k_intron=99)

        self.assertEqual(_ids(kept), {"MP_short_good"})

    def test_identity_is_weighted_by_aligned_residues(self):
        """A block aligning many residues outweighs one aligning few.

        Weighting is by query length, so a wide genomic footprint alone does
        not buy influence -- 300 residues at 0.9 against 10 at 0.1.
        """
        lines = [_seg("MP1", 100, 999, 0.9, target_span=(1, 300)),
                 _seg("MP1", 2000, 2009, 0.1, target_span=(301, 310))]

        (aln,) = summarize_alignments(lines)

        self.assertGreater(aln.identity_w, 0.85)


class TestIdentityFloor(unittest.TestCase):
    def test_floor_applies_to_alignment_not_segment(self):
        """One weak exon must not amputate an otherwise strong alignment."""
        lines = [_seg("MP1", 100, 999, 0.95), _seg("MP1", 2000, 2009, 0.20)]

        kept = cap_miniprot_lines(lines, k_base=9, k_intron=9, min_identity=0.5)

        self.assertEqual(len(kept), 2)

    def test_alignment_below_floor_is_dropped_whole(self):
        lines = [_seg("MP1", 100, 200, 0.30), _seg("MP2", 5000, 5100, 0.90)]

        kept = cap_miniprot_lines(lines, k_base=9, k_intron=9, min_identity=0.5)

        self.assertEqual(_ids(kept), {"MP2"})


class TestDeterminism(unittest.TestCase):
    def test_output_independent_of_input_order(self):
        """Input is not coordinate-sorted; results must not depend on order."""
        lines = [_seg("MP1", 100, 200, 0.9), _seg("MP2", 150, 250, 0.9),
                 _seg("MP3", 180, 280, 0.9), _seg("MP4", 900, 950, 0.7)]

        forward = cap_miniprot_lines(list(lines), k_base=2, k_intron=99)
        reverse = cap_miniprot_lines(list(reversed(lines)), k_base=2,
                                     k_intron=99)

        self.assertEqual(_ids(forward), _ids(reverse))

    def test_equal_quality_tie_broken_deterministically(self):
        lines = [_seg("MP_b", 100, 200, 0.9), _seg("MP_a", 100, 200, 0.9)]

        first = cap_miniprot_lines(list(lines), k_base=1, k_intron=99)
        second = cap_miniprot_lines(list(reversed(lines)), k_base=1,
                                    k_intron=99)

        self.assertEqual(_ids(first), _ids(second))


class TestPassthrough(unittest.TestCase):
    def test_comments_and_headers_survive(self):
        lines = ["##gff-version 3", "# note", _seg("MP1", 100, 200, 0.9)]

        kept = cap_miniprot_lines(lines, k_base=9, k_intron=9)

        self.assertEqual(kept[0], "##gff-version 3")
        self.assertEqual(kept[1], "# note")

    def test_selection_returns_ids_only(self):
        lines = [_seg("MP1", 100, 200, 0.9), _seg("MP2", 100, 200, 0.5)]

        accepted = select_alignments(summarize_alignments(lines), k_base=1,
                                     k_intron=99)

        self.assertEqual(accepted, {"MP1"})


class TestSelfOverlapCap(unittest.TestCase):
    """An alignment's own segments may overlap; its self-contribution counts.

    Checking only the *current* depth before admission lets a single ID push a
    base past the cap on its own. Real Ath miniprot output has no overlapping
    segments within an ID, so this never fired in production data -- but the
    invariant must hold for arbitrary input, not just the sample we measured.
    """

    def test_self_overlap_cannot_exceed_k_base(self):
        lines = [_seg("MP1", 100, 200, 0.9), _seg("MP1", 150, 250, 0.9)]

        kept = cap_miniprot_lines(lines, k_base=1, k_intron=99)

        self.assertLessEqual(_depth_at(kept, "C1", "+", 175), 1)

    def test_self_overlap_admitted_when_capacity_allows(self):
        lines = [_seg("MP1", 100, 200, 0.9), _seg("MP1", 150, 250, 0.9)]

        kept = cap_miniprot_lines(lines, k_base=2, k_intron=99)

        self.assertEqual(len(kept), 2)

    def test_duplicate_identical_segments_count_twice(self):
        lines = [_seg("MP1", 100, 200, 0.9), _seg("MP1", 100, 200, 0.9)]

        kept = cap_miniprot_lines(lines, k_base=1, k_intron=99)

        self.assertLessEqual(_depth_at(kept, "C1", "+", 150), 1)


class TestTargetCoordinates(unittest.TestCase):
    """Ranking must use query residues, not genomic span.

    `aligned_aa` summed genomic segment lengths, so an alignment with long
    introns outranked a compact one that aligned more of the query. The GFF
    carries the query interval (`Target=<id> <start> <end>`); use it.
    """

    def test_aligned_aa_counts_query_residues(self):
        lines = [_seg("MP1", 100, 200, 0.9, target_span=(1, 50)),
                 _seg("MP1", 4000, 4100, 0.9, target_span=(51, 100))]

        (aln,) = summarize_alignments(lines)

        self.assertEqual(aln.aligned_aa, 100)

    def test_overlapping_target_intervals_are_unioned(self):
        lines = [_seg("MP1", 100, 200, 0.9, target_span=(1, 60)),
                 _seg("MP1", 400, 500, 0.9, target_span=(41, 100))]

        (aln,) = summarize_alignments(lines)

        self.assertEqual(aln.aligned_aa, 100)

    def test_identity_weighted_by_query_length(self):
        """A segment aligning many residues dominates one aligning few, even
        when its genomic footprint is the smaller of the two."""
        lines = [_seg("MP1", 100, 130, 1.0, target_span=(1, 90)),
                 _seg("MP1", 1000, 1900, 0.0, target_span=(91, 100))]

        (aln,) = summarize_alignments(lines)

        self.assertAlmostEqual(aln.identity_w, 0.9, places=6)

    def test_reversed_query_span_is_normalised(self):
        lines = [_seg("MP1", 100, 200, 0.9, target_span=(80, 20))]

        (aln,) = summarize_alignments(lines)

        self.assertEqual(aln.aligned_aa, 61)


class TestRejectAudit(unittest.TestCase):
    """Reject reasons let us split the coverage loss by cause.

    Raw union coverage fell 28% at k=4 and only recovered 5pp at k=32, which
    points at atomicity rather than cap strength -- but that needs measuring,
    not asserting.
    """

    def test_base_cap_rejection_is_labelled(self):
        lines = [_seg("MP1", 100, 200, 0.9), _seg("MP2", 100, 200, 0.8)]

        report = audit_selection(summarize_alignments(lines), k_base=1,
                                 k_intron=99)

        self.assertEqual(report.reasons["MP2"], "base_cap")

    def test_intron_cap_rejection_is_labelled(self):
        blocks = [(100, 200), (300, 400)]
        lines = _spliced("MP1", blocks, identity=0.9) + \
            _spliced("MP2", blocks, identity=0.8)

        report = audit_selection(summarize_alignments(lines), k_base=99,
                                 k_intron=1)

        self.assertEqual(report.reasons["MP2"], "intron_cap")

    def test_identity_floor_rejection_is_labelled(self):
        lines = [_seg("MP1", 100, 200, 0.2)]

        report = audit_selection(summarize_alignments(lines), k_base=9,
                                 k_intron=9, min_identity=0.5)

        self.assertEqual(report.reasons["MP1"], "identity_floor")

    def test_accepted_alignments_have_no_reason(self):
        lines = [_seg("MP1", 100, 200, 0.9)]

        report = audit_selection(summarize_alignments(lines), k_base=1,
                                 k_intron=9)

        self.assertEqual(report.accepted, {"MP1"})
        self.assertNotIn("MP1", report.reasons)

    def test_cap_lost_coverage_counts_bases_no_longer_covered(self):
        """MP2 loses its private exon because its *other* exon hit a full base.

        MP1 and MP2 both cover 100-200; with k_base=1 only MP1 fits. MP2 also
        covers 5000-5100, which nothing else touches -- those bases end at
        depth 0 purely because of the collision elsewhere.
        """
        lines = [_seg("MP1", 100, 200, 0.9),
                 _seg("MP2", 100, 200, 0.8),
                 _seg("MP2", 5000, 5100, 0.8)]

        report = audit_selection(summarize_alignments(lines), k_base=1,
                                 k_intron=99)

        self.assertEqual(report.cap_lost_coverage_bases, 101)

    def test_no_lost_coverage_when_nothing_is_rejected(self):
        lines = [_seg("MP1", 100, 200, 0.9), _seg("MP2", 5000, 5100, 0.8)]

        report = audit_selection(summarize_alignments(lines), k_base=1,
                                 k_intron=99)

        self.assertEqual(report.cap_lost_coverage_bases, 0)


class TestAtomicityUpperBound(unittest.TestCase):
    """How much coverage could partial (non-atomic) admission possibly add?

    The conflict-ratio statistic is per alignment, so it cannot rule out "a few
    long alignments with small conflicts". This bounds the question in bp: the
    unsaturated part of every rejected alignment, unioned, minus what accepted
    alignments already cover. If that number is small, relaxing ID atomicity
    cannot recover meaningful coverage no matter how it is implemented.
    """

    def test_bound_counts_only_the_unsaturated_part(self):
        """MP2's 100-200 is full; only 201-300 could ever be rescued."""
        lines = [_seg("MP1", 100, 200, 0.9), _seg("MP2", 100, 300, 0.8)]

        report = audit_selection(summarize_alignments(lines), k_base=1,
                                 k_intron=99)

        self.assertEqual(report.rescuable_new_coverage_bases, 100)

    def test_bound_unions_overlapping_rejects(self):
        """Two rejects sharing a free region must not be counted twice."""
        lines = [_seg("MP1", 100, 200, 0.9),
                 _seg("MP2", 100, 300, 0.8),
                 _seg("MP3", 100, 300, 0.7)]

        report = audit_selection(summarize_alignments(lines), k_base=1,
                                 k_intron=99)

        self.assertEqual(report.rescuable_new_coverage_bases, 100)
        self.assertEqual(report.rescuable_bp_sum, 200)

    def test_bound_excludes_area_already_covered_by_accepted(self):
        """Unsaturated is not the same as uncovered: depth 1 < k=2 still has
        an accepted alignment on it, so rescuing there adds nothing."""
        lines = [_seg("MP1", 100, 300, 0.9),
                 _seg("MP2", 100, 200, 0.8),
                 _seg("MP3", 100, 200, 0.7)]

        report = audit_selection(summarize_alignments(lines), k_base=2,
                                 k_intron=99)

        self.assertEqual(report.rescuable_new_coverage_bases, 0)

    def test_no_bound_when_nothing_rejected(self):
        lines = [_seg("MP1", 100, 200, 0.9), _seg("MP2", 5000, 5100, 0.8)]

        report = audit_selection(summarize_alignments(lines), k_base=1,
                                 k_intron=99)

        self.assertEqual(report.rescuable_new_coverage_bases, 0)


class TestIdConsistency(unittest.TestCase):
    """An ID is one genomic hit. Anything else means the file is not what we
    think it is, and silently keeping the first value would cap the wrong
    thing -- `_target_span` drops the query name, so a mixed ID would union
    coordinates from two different proteins.
    """

    def test_same_id_on_two_seqids_is_rejected(self):
        lines = [_seg("MP1", 100, 200, 0.9, seqid="C1"),
                 _seg("MP1", 100, 200, 0.9, seqid="C2")]

        with self.assertRaises(ValueError):
            summarize_alignments(lines)

    def test_same_id_on_two_strands_is_rejected(self):
        lines = [_seg("MP1", 100, 200, 0.9, strand="+"),
                 _seg("MP1", 300, 400, 0.9, strand="-")]

        with self.assertRaises(ValueError):
            summarize_alignments(lines)

    def test_same_id_with_two_queries_is_rejected(self):
        lines = [_seg("MP1", 100, 200, 0.9, target="Q1"),
                 _seg("MP1", 300, 400, 0.9, target="Q2")]

        with self.assertRaises(ValueError):
            summarize_alignments(lines)

    def test_same_id_with_conflicting_rank_is_rejected(self):
        lines = [_seg("MP1", 100, 200, 0.9, rank=1),
                 _seg("MP1", 300, 400, 0.9, rank=2)]

        with self.assertRaises(ValueError):
            summarize_alignments(lines)


class TestMalformedInput(unittest.TestCase):
    """A row we cannot parse must not slip past the cap into EVM."""

    def test_unparseable_coordinate_is_rejected(self):
        lines = ["C1\tminiprot\tminiprot\tNOPE\t200\t100\t+\t0\tID=MP1"]

        with self.assertRaises(ValueError):
            summarize_alignments(lines)

    def test_missing_id_is_rejected(self):
        lines = ["C1\tminiprot\tminiprot\t100\t200\t100\t+\t0\tRank=1"]

        with self.assertRaises(ValueError):
            summarize_alignments(lines)

    def test_short_row_is_rejected(self):
        lines = ["C1\tminiprot\tminiprot\t100\t200"]

        with self.assertRaises(ValueError):
            summarize_alignments(lines)

    def test_non_finite_identity_is_rejected(self):
        lines = [_seg("MP1", 100, 200, float("nan"))]

        with self.assertRaises(ValueError):
            summarize_alignments(lines)

    def test_comments_still_pass_through(self):
        lines = ["##gff-version 3", "# note", _seg("MP1", 100, 200, 0.9)]

        kept = cap_miniprot_lines(lines, k_base=9, k_intron=9)

        self.assertEqual(kept[0], "##gff-version 3")


class TestPartialTargetSpan(unittest.TestCase):
    """Mixing query-coordinate and genomic-length weighting inside one
    alignment silently drops the rows that lack a Target."""

    def test_partial_target_is_rejected(self):
        good = _seg("MP1", 100, 200, 0.9, target_span=(1, 50))
        no_target = ("C1\tminiprot\tminiprot\t400\t500\t100\t+\t0\t"
                     "ID=MP1;Rank=1;Identity=0.9")

        with self.assertRaises(ValueError):
            summarize_alignments([good, no_target])

    def test_uniformly_absent_target_falls_back_to_genomic(self):
        rows = [("C1\tminiprot\tminiprot\t100\t200\t100\t+\t0\t"
                 "ID=MP1;Rank=1;Identity=0.9")]

        (aln,) = summarize_alignments(rows)

        self.assertEqual(aln.aligned_aa, 101)


class TestCapValidation(unittest.TestCase):
    def test_zero_k_base_is_rejected(self):
        with self.assertRaises(ValueError):
            select_alignments([], k_base=0, k_intron=1)

    def test_negative_k_intron_is_rejected(self):
        with self.assertRaises(ValueError):
            select_alignments([], k_base=1, k_intron=-1)


class TestFileInterface(unittest.TestCase):
    """The production path streams through external sort; it must agree with
    the in-memory path and must not leave a partial file behind on failure."""

    def _write(self, path, lines):
        with open(path, "w") as fh:
            for line in lines:
                fh.write(line + "\n")

    def test_file_path_matches_in_memory_result(self):
        import tempfile
        lines = [_seg("MP1", 100, 200, 0.9), _seg("MP2", 150, 250, 0.8),
                 _seg("MP3", 180, 280, 0.7), _seg("MP4", 900, 950, 0.95)]
        expected = _ids(cap_miniprot_lines(lines, k_base=2, k_intron=9))

        with tempfile.TemporaryDirectory() as tmp:
            src, dst = f"{tmp}/in.gff", f"{tmp}/out.gff"
            self._write(src, lines)
            cap_miniprot_file(src, dst, k_base=2, k_intron=9)
            with open(dst) as fh:
                got = _ids([ln.rstrip("\n") for ln in fh if ln.strip()])

        self.assertEqual(got, expected)

    def test_unsorted_input_is_handled(self):
        """Real files are not coordinate sorted and interleave IDs."""
        import tempfile
        lines = [_seg("MP1", 100, 200, 0.9), _seg("MP2", 5000, 5100, 0.8),
                 _seg("MP1", 400, 500, 0.9), _seg("MP2", 5400, 5500, 0.8)]

        with tempfile.TemporaryDirectory() as tmp:
            src, dst = f"{tmp}/in.gff", f"{tmp}/out.gff"
            self._write(src, lines)
            cap_miniprot_file(src, dst, k_base=9, k_intron=9)
            with open(dst) as fh:
                got = [ln for ln in fh if ln.strip()]

        self.assertEqual(len(got), 4)

    def test_failure_leaves_previous_output_intact(self):
        import tempfile
        bad = ["C1\tminiprot\tminiprot\tNOPE\t200\t100\t+\t0\tID=MP1"]

        with tempfile.TemporaryDirectory() as tmp:
            src, dst = f"{tmp}/in.gff", f"{tmp}/out.gff"
            self._write(src, bad)
            self._write(dst, ["PREVIOUS"])
            with self.assertRaises(ValueError):
                cap_miniprot_file(src, dst, k_base=1, k_intron=1)
            with open(dst) as fh:
                self.assertEqual(fh.read().strip(), "PREVIOUS")


if __name__ == "__main__":
    unittest.main()
