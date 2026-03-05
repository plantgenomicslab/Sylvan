#!/usr/bin/env bash

set -x
set -e

mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"
export SLURM_TMPDIR="$TMPDIR"

# Pipeline config: input paths, cutoffs, thread counts (read by Snakefiles)
export SYLVAN_FILTER_CONFIG="toydata/config/config_filter.yml"

# Cluster config: SLURM account, partition, per-rule resources (read by Snakemake --cluster-config)
export SYLVAN_FILTER_CLUSTER_CONFIG="toydata/config/cluster_filter.yml"

# Cluster submit command — account and partition are optional (skipped when empty or 'placeholder')
CLUSTER_CMD="python3 bin/cluster_submit.py {cluster.nodes} {cluster.memory} {cluster.ncpus} {cluster.name} {cluster.account} {cluster.partition} {cluster.time} {cluster.output} {cluster.error} {cluster.extra_args}"

# Print log location on exit (success or failure)
trap 'echo ""; echo "=== Log files: results/logs/{rule}_{wildcards}.err ==="; echo "Debug: cat results/logs/RULENAME_*.err"; echo "Recent: ls -lt results/logs/*.err | head"' EXIT

snakemake -p \
	--rerun-incomplete \
	--use-singularity \
	--singularity-args "--nv -B /data/gpfs" \
	--keep-going \
	--keep-incomplete \
	--stats filter_runtime_stats.json \
	--cluster-config "$SYLVAN_FILTER_CLUSTER_CONFIG" \
	--snakefile bin/Snakefile_filter \
	--max-jobs-per-second 50 \
	--max-status-checks-per-second 50 \
	--jobs 150 \
	--latency-wait 30 \
	--cluster "$CLUSTER_CMD" \
		"$@"

# Generate report after run
SYLVAN_FILTER_CONFIG="toydata/config/config_filter.yml" snakemake --report results/report_filter.html --snakefile bin/Snakefile_filter

# To force rerun:
# ./bin/filter_toydata.sh --forceall
