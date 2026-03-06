#!/usr/bin/env bash
#
# Local filter pipeline runner (no SLURM)
# Runs Snakemake with --cores instead of --cluster
#
# Usage:
#   ./bin/filter_local.sh              # normal run
#   ./bin/filter_local.sh -n           # dry-run
#   ./bin/filter_local.sh --forceall   # force rerun all rules
#

set -x
set -e

mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"

# Pipeline config: input paths, cutoffs, thread counts
export SYLVAN_FILTER_CONFIG="${SYLVAN_FILTER_CONFIG:-toydata/config/config_filter_local.yml}"

# Singularity args: --nv enables NVIDIA GPU passthrough (safe to include even without GPU)
REAL_PWD="$(pwd -P)"
SINGULARITY_ARGS="${SYLVAN_SINGULARITY_ARGS:---nv -B $(pwd) -B ${REAL_PWD} -B /tmp}"
export SNAKEMAKE_SINGULARITY_ARGS="$SINGULARITY_ARGS"

# Verify config file exists
if [ ! -f "$SYLVAN_FILTER_CONFIG" ]; then
	echo "ERROR: Config file not found: $SYLVAN_FILTER_CONFIG" >&2
	exit 1
fi

# Print log location on exit (success or failure)
trap 'echo ""; echo "=== Log files: results/logs/{rule}_{wildcards}.err ==="; echo "Debug: cat results/logs/RULENAME_*.err"; echo "Recent: ls -lt results/logs/*.err | head"' EXIT

# Local execution: no --cluster, use --cores for local parallelism
snakemake -p \
	--rerun-incomplete \
	--use-singularity \
	--singularity-args "$SINGULARITY_ARGS" \
	--keep-going \
	--keep-incomplete \
	--stats filter_runtime_stats.json \
	--snakefile bin/Snakefile_filter \
	--cores 16 \
	--latency-wait 10 \
		"$@"
