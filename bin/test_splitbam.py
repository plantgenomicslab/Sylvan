"""Unit tests for splitBam bin planning (psiClass per-chromosome + dense-region splitting, #25).

compute_bins packs sequences into coarse bins (one region per whole sequence). subdivide_bins
then splits any single-chromosome bin wider than max_span into sub-region bins so no region hands
psiClass a subexon graph big enough to overflow its int32 atCnt counter and hang. plan_subregions
is the pure tiling core: it cuts at zero-coverage gaps where available and at a forced boundary
otherwise (the only option on a gaplessly covered arm like A. thaliana Chr1).
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# pylint: disable=wrong-import-position
from splitBam import (  # noqa: E402
    PROFILE_BIN,
    compute_bins,
    make_forced_cut_selector,
    plan_subregions,
    subdivide_bins,
)


# ---- compute_bins: coarse binning, now region-shaped (binid, [(seq, start, end)]) ----

def test_no_big_single_bin():
    assert compute_bins([("s1", 1000), ("s2", 2000)], threshold=5_000_000) == \
        [("all", [("s1", 1, 1000), ("s2", 1, 2000)])]


def test_big_each_own_bin_full_length_region():
    b = compute_bins([("chr1", 10_000_000), ("chr2", 8_000_000)], threshold=5_000_000)
    assert b == [("chr1", [("chr1", 1, 10_000_000)]), ("chr2", [("chr2", 1, 8_000_000)])]


def test_mixed_small_packed_to_min_big():
    seqs = [("chr1", 10_000_000), ("chr2", 6_000_000),
            ("s1", 3_000_000), ("s2", 3_000_000), ("s3", 2_000_000)]
    b = compute_bins(seqs, threshold=5_000_000)
    sizes = dict(seqs)
    small_bins = [x for x in b if x[0].startswith("small_")]
    assert len([x for x in b if x[0] in ("chr1", "chr2")]) == 2   # big individual
    # each small bin's summed sequence length <= smallest big chromosome
    assert all(sum(sizes[seq] for seq, _, _ in regs) <= 6_000_000 for _, regs in small_bins)
    assert sorted(len(regs) for _, regs in small_bins) == [1, 2]  # 3+3 packed, 2 alone


def test_no_sequence_lost():
    seqs = [("chr1", 20_000_000)] + [(f"s{i}", 1_000_000) for i in range(50)]
    b = compute_bins(seqs, threshold=5_000_000)
    covered = [seq for _, regs in b for seq, _, _ in regs]
    assert sorted(covered) == sorted(n for n, _ in seqs)


def test_sanitize_binid_keeps_original_region_name():
    b = compute_bins([("chr|1", 10_000_000)], threshold=5_000_000)
    assert b == [("chr_1", [("chr|1", 1, 10_000_000)])]  # filename-safe id, original seq name


def test_sanitized_name_collision_disambiguated():
    b = compute_bins([("chr:1", 10_000_000), ("chr|1", 8_000_000)], threshold=5_000_000)
    assert [bid for bid, _ in b] == ["chr_1", "chr_1_1"]


# ---- plan_subregions: pure tiling of [1, length] into <= max_span sub-regions ----

def test_plan_single_region_when_within_budget():
    assert plan_subregions(4_000_000, [], 5_000_000) == [(1, 4_000_000)]
    assert plan_subregions(5_000_000, [], 5_000_000) == [(1, 5_000_000)]


def test_plan_gapless_arm_tiles_bounded_and_lossless():
    # A. thaliana Chr1's hung first arm is 12.8Mb of continuous coverage (zero gaps).
    regs = plan_subregions(12_827_620, [], 5_000_000)
    assert len(regs) == 3                       # ceil(12.83 / 5)
    assert regs[0][0] == 1 and regs[-1][1] == 12_827_620
    for start, end in regs:
        assert 1 <= end - start + 1 <= 5_000_000
    for (_, e0), (s1, _) in zip(regs, regs[1:]):
        assert s1 == e0 + 1                     # tiles with no base lost, no overlap


def test_plan_cuts_inside_gap_when_available():
    # a zero-coverage gap just under the 5Mb boundary -> the cut lands inside the gap
    regs = plan_subregions(9_000_000, [(4_900_000, 4_960_000)], 5_000_000)
    assert 4_900_000 <= regs[0][1] <= 4_960_000  # first sub-region ends within the gap
    for (_, e0), (s1, _) in zip(regs, regs[1:]):
        assert s1 == e0 + 1


def test_plan_prefers_gap_closest_to_boundary():
    # two gaps in reach; the one nearer the boundary is chosen to maximize the sub-region
    regs = plan_subregions(9_000_000, [(1_000_000, 1_001_000), (4_800_000, 4_802_000)], 5_000_000)
    assert regs[0][1] == 4_802_000              # cut at end of the later gap, not the early one


def test_plan_gap_spanning_boundary_cuts_at_boundary_safely():
    # a gap straddling the max_span boundary -> cut exactly at the boundary (still zero-cov)
    regs = plan_subregions(9_000_000, [(4_990_000, 5_100_000)], 5_000_000)
    assert regs[0][1] == 5_000_000


def test_plan_forced_selector_steers_cut_when_no_gap():
    calls = []

    def selector(lo, hi):
        calls.append((lo, hi))
        return lo + 1234                        # a "thin" point near the region start

    regs = plan_subregions(9_000_000, [], 5_000_000, forced_cut_selector=selector)
    assert calls                                # selector consulted (no safe gaps)
    assert regs[0] == (1, 1235)


def test_plan_no_base_lost_property_randomized():
    random.seed(0)
    for _ in range(200):
        length = random.randint(1, 30_000_000)
        max_span = random.randint(100_000, 6_000_000)
        gaps, pos = [], 1
        while pos < length:
            pos += random.randint(50_000, 3_000_000)
            glen = random.randint(1, 8_000)
            if pos + glen < length:
                gaps.append((pos, pos + glen))
                pos += glen
        regs = plan_subregions(length, gaps, max_span)
        assert regs[0][0] == 1
        assert regs[-1][1] == length            # full coverage of [1, length]
        for start, end in regs:
            assert end >= start
            assert end - start + 1 <= max_span  # bounded -> psiClass component bounded
        for (_, e0), (s1, _) in zip(regs, regs[1:]):
            assert s1 == e0 + 1                 # contiguous, non-overlapping tiling


# ---- make_forced_cut_selector: forced cuts land at the thinnest window ----

def test_forced_cut_selector_picks_min_depth_window():
    # windows of PROFILE_BIN bp; window index 3 is the thinnest in the search range
    profile = [100, 90, 80, 5, 85, 95, 100]
    selector = make_forced_cut_selector(profile, PROFILE_BIN)
    cut = selector(1, 7 * PROFILE_BIN)
    assert 3 * PROFILE_BIN <= cut <= 4 * PROFILE_BIN   # inside the min-depth window
    assert cut == 3 * PROFILE_BIN + PROFILE_BIN // 2 + 1


def test_forced_cut_selector_clamps_to_window():
    selector = make_forced_cut_selector([50, 50, 50], PROFILE_BIN)
    lo, hi = 12_345, 18_000
    cut = selector(lo, hi)
    assert lo <= cut <= hi


# ---- subdivide_bins: only over-large single-chromosome bins are split ----

def _fake_gapless_scan(_seq, length):
    """Coverage stub: no zero-cov gaps, flat profile (worst case -> all forced cuts)."""
    return [], [0] * (length // PROFILE_BIN + 1)


def test_subdivide_splits_only_oversize_chromosome_bins():
    coarse = compute_bins([("chr1", 12_827_620), ("chr2", 6_000_000)], threshold=5_000_000)
    bins = subdivide_bins("x.bam", coarse, max_span=5_000_000, scan=_fake_gapless_scan)
    chr1_bins = [(bid, regs) for bid, regs in bins if bid.startswith("chr1")]
    chr2_bins = [(bid, regs) for bid, regs in bins if bid.startswith("chr2")]
    assert [bid for bid, _ in chr1_bins] == ["chr1_r0", "chr1_r1", "chr1_r2"]
    assert [bid for bid, _ in chr2_bins] == ["chr2_r0", "chr2_r1"]
    for _, regs in chr1_bins + chr2_bins:
        (_, start, end), = regs                 # exactly one region per sub-bin
        assert end - start + 1 <= 5_000_000     # every bin bounded


def test_subdivide_is_lossless_across_subregions():
    coarse = compute_bins([("chr1", 12_827_620)], threshold=5_000_000)
    bins = subdivide_bins("x.bam", coarse, max_span=5_000_000, scan=_fake_gapless_scan)
    regs = sorted((s, e) for _, rr in bins for (_, s, e) in rr)
    assert regs[0][0] == 1
    assert regs[-1][1] == 12_827_620            # tiles the whole chromosome
    for (_, e0), (s1, _) in zip(regs, regs[1:]):
        assert s1 == e0 + 1                     # no base lost, no overlap across bins


def test_subdivide_leaves_small_and_in_budget_bins_untouched():
    # small-scaffold bins (many seqs) and chromosomes <= max_span pass through unchanged
    coarse = compute_bins(
        [("chr1", 6_000_000), ("chrSmall", 4_000_000),
         ("s1", 1_000_000), ("s2", 1_000_000)],
        threshold=5_000_000,
    )
    bins = subdivide_bins("x.bam", coarse, max_span=8_000_000, scan=_fake_gapless_scan)
    # chr1 (6Mb) <= 8Mb budget, chrSmall + s1/s2 packed into a multi-seq small bin
    assert ("chr1", [("chr1", 1, 6_000_000)]) in bins
    assert any(bid.startswith("small_") and len(regs) > 1 for bid, regs in bins)
    assert not any("_r" in bid for bid, _ in bins)   # nothing was subdivided


def test_subdivide_cuts_at_gap_when_present():
    def scan_with_gap(_seq, length):
        # one wide zero-cov gap near 5Mb; flat profile elsewhere
        return [(4_900_000, 4_970_000)], [0] * (length // PROFILE_BIN + 1)

    coarse = compute_bins([("chr1", 9_000_000)], threshold=5_000_000)
    bins = subdivide_bins("x.bam", coarse, max_span=5_000_000, scan=scan_with_gap)
    first_end = bins[0][1][0][2]
    assert 4_900_000 <= first_end <= 4_970_000       # cut fell inside the real gap


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS: {fn.__name__}")
    print(f"\n{len(fns)} tests passed")
