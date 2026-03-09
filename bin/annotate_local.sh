#!/usr/bin/env bash
#
# Local annotation pipeline runner (no SLURM)
# Runs Snakemake with --cores instead of --cluster
#
# Usage:
#   ./bin/annotate_local.sh              # normal run
#   ./bin/annotate_local.sh -n           # dry-run
#   ./bin/annotate_local.sh --forceall   # force rerun all rules
#

set -x
set -e

mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"

# Pipeline config: input paths, tool parameters, thread counts
export SYLVAN_CONFIG="${SYLVAN_CONFIG:-toydata/config/config_annotate_local.yml}"

# Singularity args: --nv enables NVIDIA GPU passthrough (safe to include even without GPU)
# Bind both symlink path and resolved real path for container access
REAL_PWD="$(pwd -P)"
SINGULARITY_ARGS="${SYLVAN_SINGULARITY_ARGS:---nv -B $(pwd) -B ${REAL_PWD} -B /tmp}"
export SNAKEMAKE_SINGULARITY_ARGS="$SINGULARITY_ARGS"

# Verify config file exists
if [ ! -f "$SYLVAN_CONFIG" ]; then
	echo "ERROR: Config file not found: $SYLVAN_CONFIG" >&2
	exit 1
fi

# Print log location on exit (success or failure)
trap 'echo ""; echo "=== Log files: results/logs/{rule}_{wildcards}.err ==="; echo "Debug: cat results/logs/RULENAME_*.err"; echo "Recent: ls -lt results/logs/*.err | head"' EXIT

# Local execution: no --cluster, use --cores for local parallelism
# Singularity bind paths adapted for local filesystem
snakemake -p \
	--rerun-incomplete \
	--rerun-triggers mtime \
	--use-singularity \
	--singularity-args "$SINGULARITY_ARGS" \
	--keep-going \
	--keep-incomplete \
	--stats annotation_runtime_stats.json \
	--snakefile bin/Snakefile_annotate \
	--groups Sam2Transfrag=group0 --group-components group0=100 \
	--cores 16 \
	--latency-wait 120 \
		"$@"
