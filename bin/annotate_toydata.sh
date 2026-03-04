#!/usr/bin/env bash
#
# Toydata annotation pipeline runner
# SLURM: cpu-s1-pgl-0 (4 nodes, 64 cores/node, 256GB/node, max 14 days)
#
# Usage:
#   ./bin/annotate_toydata.sh              # normal run
#   ./bin/annotate_toydata.sh -n           # dry-run
#   ./bin/annotate_toydata.sh --forceall   # force rerun all rules
#
# To regenerate cluster config with auto-time detection:
#   python3 bin/generate_cluster_from_config.py \
#     --config toydata/config/config_annotate.yml \
#     --out toydata/config/cluster_annotate.yml \
#     --account cpu-s1-pgl-0 --partition cpu-s1-pgl-0

set -x
set -e

mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"
export SLURM_TMPDIR="$TMPDIR"

# Pipeline config: input paths, tool parameters, thread counts (read by Snakefiles)
export SYLVAN_CONFIG="toydata/config/config_annotate.yml"

# Cluster config: SLURM account, partition, per-rule resources (read by Snakemake --cluster-config)
export SYLVAN_CLUSTER_CONFIG="toydata/config/cluster_annotate.yml"

# Verify config files exist
for cfg in "$SYLVAN_CONFIG" "$SYLVAN_CLUSTER_CONFIG"; do
	if [ ! -f "$cfg" ]; then
		echo "ERROR: Config file not found: $cfg" >&2
		exit 1
	fi
done

# Cluster submit command — account and partition are optional (skipped when empty or 'placeholder')
CLUSTER_CMD="python3 bin/cluster_submit.py {cluster.nodes} {cluster.memory} {cluster.ncpus} {cluster.name} {cluster.account} {cluster.partition} {cluster.time} {cluster.output} {cluster.error} {cluster.extra_args}"

# Print log location on exit (success or failure)
trap 'echo ""; echo "=== Log files: results/logs/{rule}_{wildcards}.err ==="; echo "Debug: cat results/logs/RULENAME_*.err"; echo "Recent: ls -lt results/logs/*.err | head"' EXIT

snakemake -p \
	--rerun-incomplete \
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

# To generate report after run:
# snakemake --report results/report.html --snakefile bin/Snakefile_annotate

# To force rerun:
# ./bin/annotate_toydata.sh --forceall
