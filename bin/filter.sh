#!/usr/bin/env bash

set -x
set -e

mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"
export SLURM_TMPDIR="$TMPDIR"

# Set config path (can be overridden via SYLVAN_FILTER_CONFIG environment variable)
export SYLVAN_FILTER_CONFIG="${SYLVAN_FILTER_CONFIG:-config/config_filter.yml}"

# Verify config file exists
if [ ! -f "$SYLVAN_FILTER_CONFIG" ]; then
	echo "ERROR: Config file not found: $SYLVAN_FILTER_CONFIG" >&2
	exit 1
fi

# Print log location on exit (success or failure)
trap 'echo ""; echo "=== Log files: results/logs/{rule}_{wildcards}.err ==="; echo "Debug: cat results/logs/RULENAME_*.err"; echo "Recent: ls -lt results/logs/*.err | head"' EXIT

snakemake -p --rerun-incomplete --cluster-config "$SYLVAN_FILTER_CONFIG" \
		--rerun-triggers mtime \
		--snakefile bin/Snakefile_filter \
		--use-singularity \
		--use-conda \
		--keep-going \
		--keep-incomplete \
		--stats filter_runtime_stats.json \
		--max-jobs-per-second 50 \
		--max-status-checks-per-second 50 \
		--jobs 150 \
		--latency-wait 30 \
		--cluster "sbatch -N {cluster.nodes} --mem={cluster.memory} --cpus-per-task={cluster.ncpus} \
				-J {cluster.name} \
				--parsable -A {cluster.account} -p {cluster.partition} \
				-t {cluster.time} -o {cluster.output} -e {cluster.error}" \
		--singularity-args "--cleanenv --env PYTHONNOUSERSITE=1" \
		"$@"
