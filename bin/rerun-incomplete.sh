#!/usr/bin/env bash

set -x
set -e

# Set config path (can be overridden via SYLVAN_CONFIG environment variable)
export SYLVAN_CONFIG="${SYLVAN_CONFIG:-config/config_annotate.yml}"
# Derive cluster config from SYLVAN_CONFIG
CLUSTER_CONFIG="${SYLVAN_CLUSTER_CONFIG:-$(dirname "$SYLVAN_CONFIG")/cluster_annotate.yml}"

# Rerun incomplete jobs with mtime trigger to detect changed files
snakemake -p \
	--rerun-incomplete \
	--rerun-triggers mtime \
	--use-singularity \
	--use-conda \
	--keep-going \
	--keep-incomplete \
	--stats annotation_runtime_stats.json \
	--cluster-config "$CLUSTER_CONFIG" \
	--snakefile bin/Snakefile_annotate \
	--groups Sam2Transfrag=group0 --group-components group0=100 \
	--max-jobs-per-second 50 \
	--max-status-checks-per-second 50 \
	--jobs 150 \
	--latency-wait 30 \
	--cluster "sbatch --mem={cluster.memory} --cpus-per-task={cluster.ncpus} \
			-J {cluster.name} \
			--parsable -A {cluster.account} -p {cluster.partition} \
			-t {cluster.time} -o {cluster.output} -e {cluster.error}" \
	--singularity-args "--cleanenv --env PYTHONNOUSERSITE=1" \
		"$@"

# To generate report after run:
# snakemake --report report.html --snakefile bin/Snakefile_annotate
