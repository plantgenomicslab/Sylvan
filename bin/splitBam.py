import os
import re
import subprocess
import sys

# Chromosomes longer than this get their own psiClass bin; shorter scaffolds are
# greedily bin-packed into bins whose total size is ~the smallest big chromosome,
# so every parallel psiClass job carries a comparable load (#25).
THRESHOLD = 5_000_000

# psiClass's `classes` assembler enumerates every candidate transcript path through
# the subexon graph of each connected component ("gene") and counts them in an int32
# (`atCnt`). On a whole chromosome arm the component reaches 40k+ subexons, `atCnt`
# overflows 2^31 (observed atCnt=-1601035960 on A. thaliana Chr1, seCnt=46679) and the
# subsequent DP over that graph never terminates (18h+ single-threaded). Per-chromosome
# binning (the original #25 fix) is not enough: Chr1/Chr2/Chr3 each hold ONE component
# spanning a whole 12-14Mb arm. psiClass has no cap on component size (for seCnt>10000 it
# just swaps its O(seCnt^2) DP table for a hash table and keeps going), and
# --maxDpConstraintSize truncates constraints without shrinking the component, so it
# cannot bound the worst case. The only robust bound is to keep each region fed to
# psiClass small enough that its largest component stays well under the tractable
# ceiling. A. thaliana Chr4 completed at seCnt=31445 over 12.6Mb (~2500 subexons/Mb);
# Chr3 hung at seCnt=41200 over only 8.2Mb (~5000/Mb). Capping the covered span at ~5Mb
# keeps even the densest observed regions under ~25k subexons (below the 31k that
# completed), with margin. Configurable via config psiClass_bin.max_region_span.
MAX_REGION_SPAN = 5_000_000

# Coverage scan settings for sub-region splitting (see split logic below).
MAPQ = 10             # match psiClass_bin's `samtools view -bq 10` filter
MIN_GAP = 50          # a zero-coverage run >= this many bp is a candidate (safe) cut zone
PROFILE_BIN = 10_000  # window size for the min-depth profile used to steer forced cuts
REFINE_WINDOW = 250_000  # search this far below the max_span boundary for a thin forced cut


def check_command(command):
    try:
        subprocess.run([command, '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        print(f"Error: {command} is not installed or not working properly.")
        sys.exit(1)


def check_and_index_bai(input_bam):
    bai_file = input_bam + '.bai'
    stale = (not os.path.exists(bai_file) or os.path.getsize(bai_file) == 0
             or os.path.getmtime(bai_file) < os.path.getmtime(input_bam))
    if stale:
        print(f"Index file {bai_file} is missing, empty, or outdated. Reindexing...")
        subprocess.run(['sambamba', 'index', '--nthreads=4', input_bam], check=True)
    else:
        print(f"Index file {bai_file} is up-to-date.")


def get_sequences(input_bam):
    """Return [(name, length), ...] from the BAM @SQ header lines, in header order."""
    result = subprocess.run(['samtools', 'view', '-H', input_bam], check=True, stdout=subprocess.PIPE)
    lines = result.stdout.decode().splitlines()
    sequences = []
    for line in lines:
        if not line.startswith('@SQ'):
            continue
        name = length = None
        for part in line.split('\t'):
            if part.startswith('SN:'):
                name = part[3:]
            elif part.startswith('LN:'):
                length = int(part[3:])
        if name is not None and length is not None:
            sequences.append((name, length))
    return sequences


def sanitize(name):
    """Filename-safe bin id derived from a sequence name (regions keep original names)."""
    return re.sub(r'[^A-Za-z0-9._-]', '_', name)


def compute_bins(sequences, threshold=THRESHOLD):
    """Group sequences into coarse psiClass bins (one region per whole sequence).

    sequences: list of (name, length).
    Returns: list of (binid, [(seq, start, end)]) -- pure, no I/O. Each region here is
    the full sequence (seq, 1, length); over-large per-chromosome bins are split into
    smaller sub-region bins later by subdivide_bins() using coverage information.

    - Each chromosome > threshold gets its own bin.
    - Remaining (<= threshold) scaffolds are greedily bin-packed into bins whose
      total length <= the smallest big chromosome (load balancing). Such bins hold
      many independent sequences, each <= threshold, so their largest connected
      component is already bounded -- only big single-chromosome bins can explode.
    - If no chromosome exceeds threshold, everything goes in one bin ("all"),
      i.e. no parallelization for small/fragmented genomes.
    """
    big = [(n, l) for n, l in sequences if l > threshold]
    small = [(n, l) for n, l in sequences if l <= threshold]

    if not big:
        return [("all", [(n, 1, l) for n, l in sequences])]

    bins = []
    used = set()
    for n, l in big:
        bid = sanitize(n)
        base, k = bid, 1
        while bid in used:  # guard against a sanitized-name collision
            bid = f"{base}_{k}"
            k += 1
        used.add(bid)
        bins.append((bid, [(n, 1, l)]))

    target = min(l for _, l in big)
    small_bins = []  # each: [current_total, [(name, 1, length), ...]]
    for n, l in sorted(small, key=lambda x: -x[1]):
        placed = False
        for b in small_bins:
            if b[0] + l <= target:
                b[0] += l
                b[1].append((n, 1, l))
                placed = True
                break
        if not placed:
            small_bins.append([l, [(n, 1, l)]])
    for i, b in enumerate(small_bins):
        bins.append((f"small_{i}", b[1]))
    return bins


def plan_subregions(length, gaps, max_span, forced_cut_selector=None):
    """Tile [1, length] into sub-regions each spanning <= max_span bp -- pure, no I/O.

    length: sequence length; the region planned is [1, length] (1-based inclusive).
    gaps:   sorted list of (gstart, gend) 1-based inclusive zero-coverage runs that are
            candidate SAFE cut zones (no aligned read covers them). May be empty.
    max_span: maximum bp span of any emitted sub-region (must be >= 1).
    forced_cut_selector: optional callable(lo, hi) -> int returning a cut position in
            [lo, hi] to use when NO safe gap lies within reach of the max_span boundary.
            Lets the caller steer forced cuts toward a local coverage minimum (a likely
            intron/intergenic boundary). Defaults to `hi` (cut exactly at the boundary).

    Returns: list of (start, end) 1-based inclusive sub-regions that tile [1, length]
    with no overlap and no base lost (start_{i+1} == end_i + 1). Every sub-region has
    (end - start + 1) <= max_span. A cut falls inside a gap whenever one is available
    within a max_span step; otherwise a forced cut is made (the only option on a
    gaplessly covered chromosome arm).
    """
    if length <= 0:
        return []
    if max_span < 1:
        raise ValueError("max_span must be >= 1")
    if length <= max_span:
        return [(1, length)]

    gaps = sorted((gs, ge) for gs, ge in gaps if ge >= 1 and gs <= length)

    regions = []
    cursor = 1
    while length - cursor + 1 > max_span:
        boundary = cursor + max_span - 1  # last base this sub-region may include
        # Prefer a safe cut: the largest zero-coverage base <= boundary and >= cursor,
        # across all gaps overlapping [cursor, boundary]. This maximizes the sub-region
        # while staying safe and within max_span.
        best = None
        for gs, ge in gaps:
            if gs > boundary:
                break  # sorted: no later gap starts <= boundary
            if ge < cursor:
                continue
            cut = min(ge, boundary)
            if cut >= cursor and (best is None or cut > best):
                best = cut
        if best is not None:
            cut = best
        elif forced_cut_selector is not None:
            cut = forced_cut_selector(cursor, boundary)
            if cut < cursor or cut > boundary:  # keep it in range and progressing
                cut = boundary
        else:
            cut = boundary
        regions.append((cursor, cut))
        cursor = cut + 1
    regions.append((cursor, length))
    return regions


def scan_coverage(input_bam, seq, length, *, mapq=MAPQ, bin_size=PROFILE_BIN, min_gap=MIN_GAP):
    """Stream `samtools depth -Q mapq -r seq` once; return (gaps, profile) -- I/O.

    gaps: list of (gstart, gend) 1-based inclusive zero-coverage runs >= min_gap
          (positions no MAPQ>=mapq read covers). These are candidate safe cut zones.
    profile: per-bin minimum depth over bin_size-bp windows across [1, length]. Windows
             with no covered position read as depth 0 (thin), which is what we want for
             steering forced cuts.
    """
    nbins = length // bin_size + 1
    profile = [None] * nbins
    gaps = []
    prev = 0  # last covered position seen (0 => none yet)
    with subprocess.Popen(['samtools', 'depth', '-Q', str(mapq), '-r', seq, input_bam],
                          stdout=subprocess.PIPE) as proc:
        for raw in proc.stdout:
            parts = raw.split(b'\t')
            pos = int(parts[1])
            depth = int(parts[2])
            if pos - prev - 1 >= min_gap:
                gaps.append((prev + 1, pos - 1))  # prev==0 -> leading gap (1, pos-1)
            prev = pos
            b = (pos - 1) // bin_size
            if b < nbins:
                d = profile[b]
                if d is None or depth < d:
                    profile[b] = depth
    if proc.returncode:
        raise subprocess.CalledProcessError(proc.returncode, 'samtools depth')
    if length - prev >= min_gap:  # trailing gap (prev==0 -> whole sequence uncovered)
        gaps.append((prev + 1, length))
    profile = [0 if p is None else p for p in profile]
    return gaps, profile


def make_forced_cut_selector(profile, bin_size, refine_window=REFINE_WINDOW):
    """Return a forced_cut_selector that cuts at the thinnest window just below hi.

    Forced cuts are unavoidable where a chromosome arm is gaplessly covered; placing them
    at a local coverage minimum makes them land at a likely intron/intergenic boundary
    rather than mid-exon, minimizing the chance of splitting a real transcript across two
    bins. Only the refine_window bp immediately below the max_span boundary (hi) are
    searched, so sub-regions stay close to max_span (few bins); on flat coverage the tie
    goes to the window nearest hi, i.e. no shaving.
    """
    def selector(lo, hi):
        if not profile:
            return hi
        lo2 = max(lo, hi - refine_window + 1)
        i0 = max(0, (lo2 - 1) // bin_size)
        i1 = min(len(profile) - 1, (hi - 1) // bin_size)
        if i1 < i0:
            return hi
        best_i, best_d = i1, profile[i1]
        for i in range(i1, i0 - 1, -1):  # high->low: ties keep the window nearest hi
            if profile[i] < best_d:
                best_d, best_i = profile[i], i
        center = best_i * bin_size + bin_size // 2 + 1  # 1-based window center
        return min(max(center, lo), hi)
    return selector


def subdivide_bins(input_bam, bins, max_span=MAX_REGION_SPAN, scan=None):
    """Split any single-sequence bin whose region exceeds max_span into sub-region bins.

    Only per-chromosome big bins (one region spanning a whole chromosome) can exceed
    max_span; small-scaffold bins hold many sequences each <= THRESHOLD <= max_span, so
    their largest component is already bounded and they pass through untouched. Each
    over-large bin becomes bins "{binid}_r0", "{binid}_r1", ... whose regions tile the
    chromosome, cutting at zero-coverage gaps where available and at local coverage
    minima otherwise. `scan` (default scan_coverage) is injectable for testing.

    Returns: list of (binid, [(seq, start, end)]).
    """
    scan = scan or (lambda seq, length: scan_coverage(input_bam, seq, length))
    out = []
    for binid, regions in bins:
        if len(regions) == 1:
            seq, start, end = regions[0]
            if start == 1 and end - start + 1 > max_span:
                gaps, profile = scan(seq, end)
                selector = make_forced_cut_selector(profile, PROFILE_BIN)
                subs = plan_subregions(end, gaps, max_span, selector)
                for k, (s, e) in enumerate(subs):
                    out.append((f"{binid}_r{k}", [(seq, s, e)]))
                print(f"  bin {binid}: {end} bp -> {len(subs)} sub-regions "
                      f"({len(gaps)} zero-cov gaps >= {MIN_GAP}bp)")
                continue
        out.append((binid, regions))
    return out


def write_bins(input_bam, bins, outdir):
    for binid, regions in bins:
        output_bam = os.path.join(outdir, f"bin_{binid}.bam")
        region_args = [f"{seq}:{start}-{end}" for seq, start, end in regions]
        print(f"Writing bin {binid} ({len(regions)} region) -> {output_bam}")
        with open(output_bam, 'wb') as bam_fh:
            subprocess.run(['samtools', 'view', '-b', input_bam, *region_args], check=True, stdout=bam_fh)


def main():
    if len(sys.argv) not in (3, 4):
        print("Usage: python splitBam.py <input_bam> <output_dir> [max_region_span]")
        sys.exit(1)

    input_bam = sys.argv[1]
    output_dir = sys.argv[2]
    max_span = int(sys.argv[3]) if len(sys.argv) == 4 else MAX_REGION_SPAN

    check_command('samtools')
    check_command('sambamba')
    check_and_index_bai(input_bam)

    os.makedirs(output_dir, exist_ok=True)
    sequences = get_sequences(input_bam)
    coarse = compute_bins(sequences)
    bins = subdivide_bins(input_bam, coarse, max_span)
    print(f"{len(sequences)} sequences -> {len(coarse)} coarse bins -> {len(bins)} bins "
          f"(threshold {THRESHOLD} bp, max region span {max_span} bp)")
    write_bins(input_bam, bins, output_dir)


if __name__ == "__main__":
    main()
