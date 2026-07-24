#!/usr/bin/env bash

set -x
set -euo pipefail

mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"
export SLURM_TMPDIR="$TMPDIR"

# ---- fleet race guard (Sylvan-EGAPx issue #29) ---------------------------
# --rerun-incomplete below DELETES the outputs of any job Snakemake considers
# incomplete before rerunning it. In the 15-run Sylvan-EGAPx fleet each run is
# driven by a self-perpetuating controller (SLURM job "syl_an_<RUN>" for annotate
# or "syl_fi_<RUN>" for filter) that runs `snakemake --unlock` at the start of
# every wave. If a manual rerun-incomplete acquires the Snakemake lock while a
# controller is between waves, the next wave's --unlock forcibly breaks that lock
# and both Snakemake instances submit jobs against the same outputs (duplicate
# DAG / concurrent writers). Because the controller pre-submits its successor with
# `--dependency=afterany`, a PENDING or RUNNING "syl_(an|fi)_<RUN>" job sits in the
# queue throughout the run's active life (including the inter-wave gap), so
# matching it here refuses to run while a controller owns this run.
#
# Harmless no-op outside the fleet: on clusters/accounts without these job names
# (or without squeue) nothing matches and the script proceeds. Pass --force to
# override -- only when the controller is confirmed dead AND you have verified no
# orphaned child jobs are still RUNNING for this run.
_force=0
_passthru=()
for _a in "$@"; do
    if [ "$_a" = "--force" ]; then _force=1; else _passthru+=("$_a"); fi
done
set -- ${_passthru[@]+"${_passthru[@]}"}

_run="$(basename "$PWD")"
if [ "$_force" -ne 1 ] \
   && squeue -A cpu-s1-pgl-0 -h -o '%200j' 2>/dev/null | grep -qE "^syl_(an|fi)_${_run}$"; then
    echo "ABORT (#29): controller syl_an_${_run}/syl_fi_${_run} is queued or running for '${_run}'." >&2
    echo "  Its per-wave 'snakemake --unlock' races this manual rerun-incomplete and can spawn a" >&2
    echo "  duplicate DAG writing the same outputs. Wait for the controller chain to finish, or --" >&2
    echo "  if it is confirmed dead and no orphaned child jobs are RUNNING -- re-run with --force." >&2
    exit 1
fi
# --------------------------------------------------------------------------

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
	--singularity-args "--nv -B /data/gpfs" \
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
	--cluster-status "python3 bin/cluster_status.py" \
		"$@"

# Generate report after run
snakemake --report results/report.html --snakefile bin/Snakefile_annotate
