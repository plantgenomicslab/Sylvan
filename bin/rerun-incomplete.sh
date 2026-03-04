#!/usr/bin/env bash

set -x
set -e

mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"
export SLURM_TMPDIR="$TMPDIR"

# Pipeline config: input paths, tool parameters, thread counts (read by Snakefiles)
export SYLVAN_CONFIG="${SYLVAN_CONFIG:-config/config_annotate.yml}"

# Cluster config: SLURM account, partition, per-rule resources (read by Snakemake --cluster-config)
# Defaults to SYLVAN_CONFIG for single-file mode. Set to a separate cluster YAML to split concerns.
export SYLVAN_CLUSTER_CONFIG="${SYLVAN_CLUSTER_CONFIG:-$SYLVAN_CONFIG}"

# Cluster submit command — account and partition are optional (skipped when empty or 'placeholder')
CLUSTER_CMD="python3 bin/cluster_submit.py {cluster.nodes} {cluster.memory} {cluster.ncpus} {cluster.name} {cluster.account} {cluster.partition} {cluster.time} {cluster.output} {cluster.error} {cluster.extra_args}"

# Rerun incomplete jobs with mtime trigger to detect changed files
snakemake -p \
	--rerun-incomplete \
	--rerun-triggers mtime \
	--use-singularity \
	--use-conda \
	--keep-going \
	--keep-incomplete \
	--stats annotation_runtime_stats.json \
	--cluster-config "$SYLVAN_CLUSTER_CONFIG" \
	--snakefile bin/Snakefile_annotate \
	--groups Sam2Transfrag=group0 --group-components group0=100 \
	--max-jobs-per-second 50 \
	--max-status-checks-per-second 50 \
	--jobs 150 \
	--latency-wait 30 \
	--cluster "$CLUSTER_CMD" \
	--singularity-args "--cleanenv --env PYTHONNOUSERSITE=1" \
		"$@"

# Generate report after run
snakemake --report results/report.html --snakefile bin/Snakefile_annotate
