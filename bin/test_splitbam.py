"""Unit tests for splitBam.compute_bins (psiClass per-chromosome binning, #25)."""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from splitBam import compute_bins, sanitize  # noqa: E402


def test_no_big_single_bin():
    assert compute_bins([("s1", 1000), ("s2", 2000)], threshold=5_000_000) == [("all", ["s1", "s2"])]


def test_big_each_own_bin():
    b = compute_bins([("chr1", 10_000_000), ("chr2", 8_000_000)], threshold=5_000_000)
    assert b == [("chr1", ["chr1"]), ("chr2", ["chr2"])]


def test_mixed_small_packed_to_min_big():
    seqs = [("chr1", 10_000_000), ("chr2", 6_000_000),
            ("s1", 3_000_000), ("s2", 3_000_000), ("s3", 2_000_000)]
    b = compute_bins(seqs, threshold=5_000_000)
    sizes = dict(seqs)
    small_bins = [x for x in b if x[0].startswith("small_")]
    assert len([x for x in b if x[0] in ("chr1", "chr2")]) == 2   # big individual
    assert all(sum(sizes[n] for n in names) <= 6_000_000 for _, names in small_bins)  # <= min big
    assert sorted(len(n) for _, n in small_bins) == [1, 2]        # 3+3 packed, 2 alone


def test_no_sequence_lost():
    seqs = [("chr1", 20_000_000)] + [(f"s{i}", 1_000_000) for i in range(50)]
    b = compute_bins(seqs, threshold=5_000_000)
    covered = [n for _, names in b for n in names]
    assert sorted(covered) == sorted(n for n, _ in seqs)


def test_sanitize_binid_keeps_original_region():
    b = compute_bins([("chr|1", 10_000_000)], threshold=5_000_000)
    assert b == [("chr_1", ["chr|1"])]  # filename-safe id, original region name


def test_sanitized_name_collision_disambiguated():
    b = compute_bins([("chr:1", 10_000_000), ("chr|1", 8_000_000)], threshold=5_000_000)
    assert [bid for bid, _ in b] == ["chr_1", "chr_1_1"]


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS: {fn.__name__}")
    print(f"\n{len(fns)} tests passed")
