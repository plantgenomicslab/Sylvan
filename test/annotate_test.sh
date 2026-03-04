#!/usr/bin/env bash
#
# Test script: validates config parsing and Snakemake DAG using ./toydata
# Uses test/config_annotate_test.yml with reduced resource requests
#
# Usage:
#   ./test/annotate_test.sh              # dry-run (default)
#   ./test/annotate_test.sh --run        # actual SLURM submission
#   ./test/annotate_test.sh --forceall   # force rerun all rules

set -x
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_DIR"

mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"
export SLURM_TMPDIR="$TMPDIR"

# Pipeline config
export SYLVAN_CONFIG="test/config_annotate_test.yml"

# Cluster submit command
CLUSTER_CMD="python3 bin/cluster_submit.py {cluster.nodes} {cluster.memory} {cluster.ncpus} {cluster.name} {cluster.account} {cluster.partition} {cluster.time} {cluster.output} {cluster.error} {cluster.extra_args}"

# Print log location on exit
trap 'echo ""; echo "=== Log files: results/logs/{rule}_{wildcards}.err ==="; echo "Recent: ls -lt results/logs/*.err | head"' EXIT

# Check if --run flag is passed
RUN_MODE="--dry-run"
EXTRA_ARGS=()
for arg in "$@"; do
    if [ "$arg" = "--run" ]; then
        RUN_MODE=""
    else
        EXTRA_ARGS+=("$arg")
    fi
done

echo "=== Config validation ==="
python3 -c "import yaml; c=yaml.safe_load(open('$SYLVAN_CONFIG')); print('Config OK:', list(c.keys()))"

echo ""
echo "=== Snakemake ${RUN_MODE:-submission} ==="
snakemake -p \
	$RUN_MODE \
	--rerun-incomplete \
	--use-singularity \
	--use-conda \
	--keep-going \
	--keep-incomplete \
	--stats test/annotation_runtime_stats.json \
	--cluster-config "$SYLVAN_CONFIG" \
	--snakefile bin/Snakefile_annotate \
	--groups Sam2Transfrag=group0 --group-components group0=100 \
	--max-jobs-per-second 50 \
	--max-status-checks-per-second 50 \
	--jobs 150 \
	--latency-wait 30 \
	--cluster "$CLUSTER_CMD" \
	--singularity-args "--cleanenv --env PYTHONNOUSERSITE=1" \
		"${EXTRA_ARGS[@]}"
