#!/usr/bin/env bash

set -x

mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"
export SLURM_TMPDIR="$TMPDIR"

# Set config path for Snakefile
export SYLVAN_FILTER_CONFIG="toydata/config/config_filter.yml"

# Print log location on exit (success or failure)
trap 'echo ""; echo "=== Log files: results/logs/{rule}_{wildcards}.err ==="; echo "Debug: cat results/logs/RULENAME_*.err"; echo "Recent: ls -lt results/logs/*.err | head"' EXIT

snakemake -p \
	--rerun-incomplete \
	--use-singularity \
	--use-conda \
	--keep-going \
	--keep-incomplete \
	--stats filter_runtime_stats.json \
	--cluster-config toydata/config/config_filter.yml \
	--snakefile bin/Snakefile_filter \
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

# Generate report after run
SYLVAN_FILTER_CONFIG="toydata/config/config_filter.yml" snakemake --report results.html --snakefile bin/Snakefile_filter

# To force rerun:
# ./bin/filter_toydata.sh --forceall
