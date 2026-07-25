#!/usr/bin/env python3
"""Cut an oversized split-SAM into sub-regions at low-coverage points.

Background (tracker #24): sam2transfrag's transfrag search grows super-linearly
with region span. GETA's first-pass splitter cuts only where the genome has no
alignment for >=10 bp, so deep RNA-seq -- which covers the genome continuously --
leaves splits of several Mb. Measured on Osa: 103 of 1,028 splits (10%) carry
95.6% of the sequence, and every one of them was skipped by the 1 Mb guard, so
transfrag evidence was missing across almost the whole genome.

Raising sam2transfrag's own cut thresholds does not help: a 2.43 Mb split failed
to finish in 30 min at --max_expressed_base_depth 50, 200 and 1000 alike, at a
steady 31 MB RSS. It is span, not depth, that binds. Splitting that same file
into three ~810 kb pieces let all three finish (4,328 / 3,582 / 3,754 transfrags).

So: cut, and cut where it hurts least. Rather than slicing at fixed offsets --
which lands mid-gene as often as not -- each cut is placed at the lowest-coverage
bin within a window around the ideal offset, i.e. between genes wherever the data
allows. Same idea as the per-chromosome binning that fixed psiClass in #25.

Reads are assigned by POS, so every alignment lands in exactly one sub-region and
no record is duplicated (mergeTransfrag concatenates without de-duplication).

Usage: split_large_sam.py <in.sam> <out_dir> [--max-span 800000] [--bin 1000]
"""
import argparse
import os
import sys


def parse_region(path):
    """GETA names splits '<seqid>.<start>-<end>.sam'; seqids may contain dots."""
    base = os.path.basename(path)
    if base.endswith(".sam"):
        base = base[:-4]
    seqid, _, coords = base.rpartition(".")
    start, _, end = coords.partition("-")
    return seqid, int(start), int(end)


def bin_depth(sam_path, lo, hi, bin_size):
    """First pass: alignment starts per bin. Streaming, so memory is O(bins)."""
    nbins = max(1, (hi - lo) // bin_size + 1)
    depth = [0] * nbins
    n = 0
    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.split("\t", 4)
            if len(f) < 4:
                continue
            try:
                pos = int(f[3])
            except ValueError:
                continue
            idx = (pos - lo) // bin_size
            if 0 <= idx < nbins:
                depth[idx] += 1
                n += 1
    return depth, n


def choose_cuts(depth, lo, hi, bin_size, max_span, window_frac=0.25):
    """Place cuts at the emptiest bin near each ideal offset.

    Searching a window (default +/-25% of the piece length) keeps pieces close to
    the target size while still landing the cut in a coverage trough.
    """
    span = hi - lo
    n_parts = -(-span // max_span)          # ceil
    if n_parts <= 1:
        return []
    piece = span / n_parts
    half = int(piece * window_frac)
    cuts = []
    for i in range(1, n_parts):
        ideal = lo + int(piece * i)
        a, b = max(lo + 1, ideal - half), min(hi - 1, ideal + half)
        ia, ib = (a - lo) // bin_size, (b - lo) // bin_size
        if ib <= ia:
            cuts.append(ideal)
            continue
        # Tie-break on distance to the ideal offset. Without it, a window with no
        # trough (all bins equal) returns its leftmost bin and the pieces come out
        # lopsided; with it, a flat window cuts at the ideal offset instead.
        best = min(range(ia, ib + 1),
                   key=lambda k: (depth[k] if k < len(depth) else 0,
                                  abs(lo + k * bin_size - ideal)))
        cut = lo + best * bin_size + bin_size // 2
        # keep cuts strictly increasing even if two windows pick the same trough
        if cuts and cut <= cuts[-1]:
            cut = cuts[-1] + bin_size
        cuts.append(min(cut, hi - 1))
    return cuts


def write_subregions(sam_path, out_dir, seqid, bounds):
    """Second pass: route each record to the sub-region containing its POS."""
    os.makedirs(out_dir, exist_ok=True)
    paths, handles = [], []
    for i in range(len(bounds) - 1):
        p = os.path.join(out_dir, f"{seqid}.{bounds[i]}-{bounds[i+1]}.sam")
        paths.append(p)
        handles.append(open(p, "w"))
    counts = [0] * len(handles)
    try:
        with open(sam_path) as fh:
            for line in fh:
                if line.startswith("@"):
                    for h in handles:
                        h.write(line)
                    continue
                f = line.split("\t", 4)
                if len(f) < 4:
                    continue
                try:
                    pos = int(f[3])
                except ValueError:
                    continue
                # bounds is sorted and small (a handful of pieces); linear scan is fine
                for i in range(len(handles)):
                    if bounds[i] <= pos < bounds[i + 1]:
                        handles[i].write(line)
                        counts[i] += 1
                        break
                else:
                    if pos >= bounds[-1]:
                        handles[-1].write(line)
                        counts[-1] += 1
    finally:
        for h in handles:
            h.close()
    return paths, counts


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sam")
    ap.add_argument("out_dir")
    ap.add_argument("--max-span", type=int, default=800000,
                    help="target maximum sub-region span in bp (default: 800000)")
    ap.add_argument("--bin", type=int, default=1000,
                    help="coverage bin size in bp for locating cut points (default: 1000)")
    a = ap.parse_args()

    seqid, lo, hi = parse_region(a.sam)
    span = hi - lo
    if span <= a.max_span:
        print(f"{os.path.basename(a.sam)}: span {span} <= {a.max_span}, no split needed",
              file=sys.stderr)
        return 0

    depth, n_aln = bin_depth(a.sam, lo, hi, a.bin)
    cuts = choose_cuts(depth, lo, hi, a.bin, a.max_span)
    bounds = [lo] + cuts + [hi]
    paths, counts = write_subregions(a.sam, a.out_dir, seqid, bounds)

    print(f"{os.path.basename(a.sam)}: span {span/1e6:.2f} Mb, {n_aln} alignments "
          f"-> {len(paths)} sub-regions", file=sys.stderr)
    for p, c in zip(paths, counts):
        print(f"  {os.path.basename(p)}: {c} alignments", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
