#!/usr/bin/env bash

set -x
set -e

mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"
export SLURM_TMPDIR="$TMPDIR"

# Set config path (can be overridden via SYLVAN_CONFIG environment variable)
# This file serves as both pipeline config and --cluster-config for Snakemake.
export SYLVAN_CONFIG="${SYLVAN_CONFIG:-config/config_annotate.yml}"

# Rerun incomplete jobs with mtime trigger to detect changed files
snakemake -p \
	--rerun-incomplete \
	--rerun-triggers mtime \
	--use-singularity \
	--use-conda \
	--keep-going \
	--keep-incomplete \
	--stats annotation_runtime_stats.json \
	--cluster-config "$SYLVAN_CONFIG" \
	--snakefile bin/Snakefile_annotate \
	--groups Sam2Transfrag=group0 --group-components group0=100 \
	--max-jobs-per-second 50 \
	--max-status-checks-per-second 50 \
	--jobs 150 \
	--latency-wait 30 \
	--cluster "$(python3 bin/get_cluster_cmd.py "$SYLVAN_CONFIG")" \
	--singularity-args "--cleanenv --env PYTHONNOUSERSITE=1" \
		"$@"

# Generate report after run
snakemake --report results/report.html --snakefile bin/Snakefile_annotate
