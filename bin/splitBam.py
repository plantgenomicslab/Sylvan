import os
import re
import subprocess
import sys

# Chromosomes longer than this get their own psiClass bin; shorter scaffolds are
# greedily bin-packed into bins whose total size is ~the smallest big chromosome,
# so every parallel psiClass job carries a comparable load (#25).
THRESHOLD = 5_000_000


def check_command(command):
    try:
        subprocess.run([command, '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        print(f"Error: {command} is not installed or not working properly.")
        sys.exit(1)


def check_and_index_bai(input_bam):
    bai_file = input_bam + '.bai'
    if not os.path.exists(bai_file) or os.path.getsize(bai_file) == 0 or os.path.getmtime(bai_file) < os.path.getmtime(input_bam):
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
    """Group sequences into psiClass bins.

    sequences: list of (name, length).
    Returns: list of (binid, [seqnames]) -- pure, no I/O.

    - Each chromosome > threshold gets its own bin.
    - Remaining (<= threshold) scaffolds are greedily bin-packed into bins whose
      total length <= the smallest big chromosome (load balancing).
    - If no chromosome exceeds threshold, everything goes in one bin ("all"),
      i.e. no parallelization for small/fragmented genomes.
    """
    big = [(n, l) for n, l in sequences if l > threshold]
    small = [(n, l) for n, l in sequences if l <= threshold]

    if not big:
        return [("all", [n for n, _ in sequences])]

    bins = []
    used = set()
    for n, _ in big:
        bid = sanitize(n)
        base, k = bid, 1
        while bid in used:  # guard against a sanitized-name collision
            bid = f"{base}_{k}"
            k += 1
        used.add(bid)
        bins.append((bid, [n]))

    target = min(l for _, l in big)
    small_bins = []  # each: [current_total, [names]]
    for n, l in sorted(small, key=lambda x: -x[1]):
        placed = False
        for b in small_bins:
            if b[0] + l <= target:
                b[0] += l
                b[1].append(n)
                placed = True
                break
        if not placed:
            small_bins.append([l, [n]])
    for i, b in enumerate(small_bins):
        bins.append((f"small_{i}", b[1]))
    return bins


def write_bins(input_bam, bins, outdir):
    for binid, seqnames in bins:
        output_bam = os.path.join(outdir, f"bin_{binid}.bam")
        print(f"Writing bin {binid} ({len(seqnames)} seq) -> {output_bam}")
        with open(output_bam, 'wb') as bam_fh:
            subprocess.run(['samtools', 'view', '-b', input_bam, *seqnames], check=True, stdout=bam_fh)


def main():
    if len(sys.argv) != 3:
        print("Usage: python splitBam.py <input_bam> <output_dir>")
        sys.exit(1)

    input_bam = sys.argv[1]
    output_dir = sys.argv[2]

    check_command('samtools')
    check_command('sambamba')
    check_and_index_bai(input_bam)

    os.makedirs(output_dir, exist_ok=True)
    sequences = get_sequences(input_bam)
    bins = compute_bins(sequences)
    print(f"{len(sequences)} sequences -> {len(bins)} bins (threshold {THRESHOLD} bp)")
    write_bins(input_bam, bins, output_dir)


if __name__ == "__main__":
    main()
