#!/bin/bash
# Split an interleaved paired-end FASTQ into order-matched _1/_2 files.
#
# Some FASTQ are regenerated from a BAM rather than downloaded from SRA. Those
# carry BOTH mates in a single file tagged /1 and /2, shuffled relative to each
# other, and often with duplicated records (the same read emitted once per
# secondary/supplementary alignment). The filename looks single-end, so the
# paired rules never fire and RSEM rejects the transcriptome BAM with
# "has alignments with inconsistent read lengths".
#
# Approach: one disk-backed `sort` on the read name puts @X/1 immediately before
# its @X/2, so a single streaming pass can group by base id and emit the first
# /1 and first /2 of each group -- de-duplicating and mating in one go. Reads
# whose mate is missing are dropped.
#
# `seqkit pair` is deliberately not used: its own docs warn that shuffled input
# makes it buffer the unmatched remainder in RAM, and these files are 1-6 GB
# compressed.
#
# Usage: split_interleaved_fastq.sh <in.fastq.gz> <out_1.fastq.gz> <out_2.fastq.gz> [tmpdir] [sort_mem]
set -euo pipefail

IN=$1
O1=$2
O2=$3
TMP=${4:-${TMPDIR:-/var/tmp}}
SORTMEM=${5:-2G}

mkdir -p "$(dirname "$O1")" "$TMP"

# Refuse a RAM-backed tmpdir. On this cluster /tmp is a 244 GB tmpfs, so sort
# spill there counts against node memory: 20 concurrent splits filled a node to
# 235/244 GB and the OOM killer took out an unrelated pipeline controller.
fstype=$(df -T "$TMP" 2>/dev/null | awk 'NR==2{print $2}')
case "$fstype" in
	tmpfs|ramfs)
		echo "ERROR: tmpdir '$TMP' is $fstype (RAM-backed); refusing to spill sort there." >&2
		echo "       Point TMPDIR at disk-backed scratch." >&2
		exit 2;;
esac

# .part while writing so an interrupted run never leaves a truncated file that
# looks complete to Snakemake.
zcat "$IN" \
| paste - - - - \
| LC_ALL=C sort -t $'\t' -k1,1 -S "$SORTMEM" -T "$TMP" \
| awk -F'\t' -v o1="$O1.part" -v o2="$O2.part" '
	function flush() {
		if (r1 != "" && r2 != "") { print curid "\n" r1 | c1; print curid "\n" r2 | c2; npair++ }
		else if (curid != "")     { orphan++ }
	}
	BEGIN { c1 = "gzip > " o1; c2 = "gzip > " o2 }
	{
		id = $1
		mate = substr(id, length(id), 1)      # trailing 1 or 2 of "/1" | "/2"
		sub(/\/[12]$/, "", id)
		if (id != curid) { flush(); curid = id; r1 = ""; r2 = "" }
		rec = $2 "\n" $3 "\n" $4
		if      (mate == "1") { if (r1 == "") r1 = rec; else dup++ }
		else if (mate == "2") { if (r2 == "") r2 = rec; else dup++ }
		else                  { untagged++ }
	}
	END {
		flush()
		close(c1); close(c2)
		printf "pairs=%d orphans=%d duplicates_dropped=%d untagged=%d\n", npair, orphan+0, dup+0, untagged+0 > "/dev/stderr"
	}'

# Both mates must be present and equal in length, or the pair is unusable
# downstream. Checking here keeps a bad split from reaching STAR.
n1=$(zcat "$O1.part" | awk 'END{print NR/4}')
n2=$(zcat "$O2.part" | awk 'END{print NR/4}')
if [ "$n1" != "$n2" ]; then
	echo "ERROR: read count mismatch after split: _1=$n1 _2=$n2" >&2
	exit 1
fi
if [ "$n1" = "0" ]; then
	echo "ERROR: split produced zero pairs from $IN (is it really interleaved?)" >&2
	exit 1
fi

mv -f "$O1.part" "$O1"
mv -f "$O2.part" "$O2"
echo "split $(basename "$IN"): $n1 pairs" >&2
