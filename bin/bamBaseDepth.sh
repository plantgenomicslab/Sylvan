#!/usr/bin/env bash

set -e

# Temporary directory under results
mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"
export SLURM_TMPDIR="$TMPDIR"

# Calculate per-base read depth from a bam file

# Arguments:
# $1: prefix
# $2: path to bam file
# $3: threads

samtools index $2

mosdepth $1 $2 -t $3

echo "chromosome position $1" > $1.basedepth

zcat "$1.per-base.bed.gz" | \
  awk '{for(i=$2;i<$3;i++){ print $1,i,i+1,$4} }' | \
  cut -f1,3,4 -d" " >> $1.basedepth
