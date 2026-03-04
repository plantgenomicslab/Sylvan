#!/usr/bin/env python3
"""SLURM job submission wrapper for Snakemake.

Called by Snakemake via --cluster with {cluster.X} placeholders substituted.
Skips -A (account) and -p (partition) when their values are empty or 'placeholder',
so the same config works on HPC systems that don't require a SLURM account.

Any extra sbatch flags can be passed via {cluster.extra_args} in the config YAML
(e.g., extra_args: "--gres=gpu:1 --export=ALL"). The shell performs word-splitting
on extra_args before this script receives them, so multi-flag values work naturally.

Usage in entry scripts (Snakemake appends the jobscript path as the final argument):
  --cluster "python3 bin/cluster_submit.py {cluster.nodes} {cluster.memory} \
    {cluster.ncpus} {cluster.name} {cluster.account} {cluster.partition} \
    {cluster.time} {cluster.output} {cluster.error} {cluster.extra_args}"
"""

import subprocess
import sys


SKIP_VALUES = {"", "placeholder"}


def main():
    if len(sys.argv) < 11:
        print(
            "Usage: cluster_submit.py NODES MEMORY NCPUS NAME ACCOUNT PARTITION "
            "TIME OUTPUT ERROR [EXTRA_ARGS...] JOBSCRIPT",
            file=sys.stderr,
        )
        sys.exit(1)

    nodes = sys.argv[1]
    memory = sys.argv[2]
    ncpus = sys.argv[3]
    name = sys.argv[4]
    account = sys.argv[5]
    partition = sys.argv[6]
    time_limit = sys.argv[7]
    # Warn if time_limit looks like a sexagesimal artifact (large integer)
    if time_limit.isdigit() and int(time_limit) > 1000:
        print(
            f"WARNING: time_limit '{time_limit}' looks like a YAML 1.1 sexagesimal artifact. "
            f"Quote time values in your config (e.g., time: \"3-00:00:00\").",
            file=sys.stderr,
        )
    output = sys.argv[8]
    error = sys.argv[9]
    # Snakemake appends the job script as the last positional argument.
    # Everything between error and the job script is extra_args
    # (already word-split by the shell).
    jobscript = sys.argv[-1]
    extra_args = sys.argv[10:-1]

    cmd = [
        "sbatch",
        "-N", nodes,
        "--mem=" + memory,
        "--cpus-per-task=" + ncpus,
        "-J", name,
        "--parsable",
        "-t", time_limit,
        "-o", output,
        "-e", error,
    ]

    if account.lower() not in SKIP_VALUES:
        cmd.extend(["-A", account])
    if partition.lower() not in SKIP_VALUES:
        cmd.extend(["-p", partition])

    # Append any extra sbatch flags (e.g., --gres=gpu:1).
    # These are already word-split by the shell, so append directly.
    for arg in extra_args:
        if arg.strip():
            cmd.append(arg)

    cmd.append(jobscript)

    result = subprocess.run(cmd, capture_output=True, text=True)
    sys.stdout.write(result.stdout)
    if result.stderr:
        sys.stderr.write(result.stderr)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
