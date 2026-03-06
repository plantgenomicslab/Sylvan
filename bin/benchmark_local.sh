#!/usr/bin/env bash
#
# Local benchmark runner — runs BUSCO + OMArk on all configured GFF3 files
#
# Usage:
#   ./bin/benchmark_local.sh              # normal run
#   ./bin/benchmark_local.sh -n           # dry-run
#
# Requires Benchmark section in filter config with gff3_files and optional omark_db.

set -x
set -e

mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"

export SYLVAN_FILTER_CONFIG="${SYLVAN_FILTER_CONFIG:-toydata/config/config_filter_local.yml}"

# Singularity args: --nv enables NVIDIA GPU passthrough (safe to include even without GPU)
REAL_PWD="$(pwd -P)"
SINGULARITY_ARGS="${SYLVAN_SINGULARITY_ARGS:---nv -B $(pwd) -B ${REAL_PWD} -B /tmp}"
export SNAKEMAKE_SINGULARITY_ARGS="$SINGULARITY_ARGS"

if [ ! -f "$SYLVAN_FILTER_CONFIG" ]; then
	echo "ERROR: Config file not found: $SYLVAN_FILTER_CONFIG" >&2
	exit 1
fi

snakemake -p \
	--rerun-incomplete \
	--use-singularity \
	--singularity-args "$SINGULARITY_ARGS" \
	--keep-going \
	--keep-incomplete \
	--stats benchmark_runtime_stats.json \
	--snakefile bin/Snakefile_benchmark \
	--cores 16 \
	--latency-wait 10 \
		"$@"
