#!/usr/bin/env bash
#
# SLURM benchmark runner — runs BUSCO + OMArk on all configured GFF3 files
#
# Usage:
#   ./bin/benchmark.sh              # normal run
#   ./bin/benchmark.sh -n           # dry-run
#

set -x
set -e

mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"
export SLURM_TMPDIR="$TMPDIR"

export SYLVAN_FILTER_CONFIG="${SYLVAN_FILTER_CONFIG:-config/config_filter.yml}"
export SYLVAN_FILTER_CLUSTER_CONFIG="${SYLVAN_FILTER_CLUSTER_CONFIG:-$SYLVAN_FILTER_CONFIG}"

CLUSTER_CMD="python3 bin/cluster_submit.py {cluster.nodes} {cluster.memory} {cluster.ncpus} {cluster.name} {cluster.account} {cluster.partition} {cluster.time} {cluster.output} {cluster.error} {cluster.extra_args}"

for cfg in "$SYLVAN_FILTER_CONFIG" "$SYLVAN_FILTER_CLUSTER_CONFIG"; do
	if [ ! -f "$cfg" ]; then
		echo "ERROR: Config file not found: $cfg" >&2
		exit 1
	fi
done

# Singularity args: --nv enables NVIDIA GPU passthrough (safe to include even without GPU)
SINGULARITY_ARGS="${SYLVAN_SINGULARITY_ARGS:---nv -B /data/gpfs}"

trap 'echo ""; echo "=== Log files: results/logs/{rule}_{wildcards}.err ==="; echo "Recent: ls -lt results/logs/*.err | head"' EXIT

snakemake -p \
	--rerun-incomplete \
	--use-singularity \
	--singularity-args "$SINGULARITY_ARGS" \
	--keep-going \
	--keep-incomplete \
	--stats benchmark_runtime_stats.json \
	--cluster-config "$SYLVAN_FILTER_CLUSTER_CONFIG" \
	--snakefile bin/Snakefile_benchmark \
	--max-jobs-per-second 50 \
	--max-status-checks-per-second 50 \
	--jobs 50 \
	--latency-wait 30 \
	--cluster "$CLUSTER_CMD" \
		"$@"
