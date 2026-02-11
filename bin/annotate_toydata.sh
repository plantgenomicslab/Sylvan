#!/usr/bin/env bash

set -x
set -e

mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"
export SLURM_TMPDIR="$TMPDIR"

# Set config path for Snakefile
export SYLVAN_CONFIG="toydata/config/config_annotate.yml"

# Print log location on exit (success or failure)
trap 'echo ""; echo "=== Log files: results/logs/{rule}_{wildcards}.err ==="; echo "Debug: cat results/logs/RULENAME_*.err"; echo "Recent: ls -lt results/logs/*.err | head"' EXIT

#--rerun-triggers mtime \
snakemake -p \
	--rerun-incomplete \
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

# To generate report after run:
# snakemake --report results/report.html --snakefile bin/Snakefile_annotate

# To force rerun:
# ./bin/annotate_toydata.sh --forceall
