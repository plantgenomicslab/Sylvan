#!/usr/bin/env python3
"""Snakemake ``--cluster-status`` hook for SLURM.

Snakemake calls this with a single argument — the external job id that
``bin/cluster_submit.py`` printed via ``sbatch --parsable`` — and reads exactly
one of three tokens from stdout:

    running / success / failed

WHY THIS EXISTS
    Without ``--cluster-status`` Snakemake learns a child job finished ONLY from
    the ``.snakemake/tmp.<hash>/<jobid>.jobfinished|.jobfailed`` marker files the
    jobscript epilogue writes. SLURM kills that skip that epilogue — TIMEOUT,
    scancel/CANCELLED, some NODE_FAIL — leave NO marker, so the wave waits
    forever (an observed 11.9-day controller deadlock on a psiClass TIMEOUT).
    Reporting the real SLURM state turns those terminal-but-markerless jobs into
    ``failed`` → ``--keep-going`` ends the wave → the next wave retries, and the
    whole marker-deadlock class disappears.

SACCT LAG (the dangerous edge)
    A freshly submitted job may not be in the accounting DB yet, so ``sacct``
    prints nothing. Reporting ``failed`` on an empty/unknown result would kill
    live jobs, so anything we cannot positively classify as terminal is reported
    ``running`` and Snakemake re-polls. ``squeue`` (live controller state) is
    consulted first precisely to absorb this lag.

NOTES
    * Uses only ``squeue`` and ``sacct`` — both work on this cluster. (Only
      ``scontrol`` has the libreadline.so.6 problem documented for the fleet;
      this script never calls it, so no shim is needed.)
    * ``sacct -X`` restricts output to the job allocation (no ``.batch`` /
      ``.extern`` steps), giving one State line per job.
"""

import subprocess
import sys


# SLURM states that mean the job is still in flight — keep Snakemake waiting.
RUNNING_STATES = {
    "RUNNING",
    "PENDING",
    "COMPLETING",
    "CONFIGURING",
    "REQUEUED",
    "REQUEUE_HOLD",
    "REQUEUE_FED",
    "RESIZING",
    "SUSPENDED",
    "SIGNALING",
    "STAGE_OUT",
    "STOPPED",
}


def _run(cmd):
    """Run a command, returning the CompletedProcess or None on any failure."""
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None


def classify(jobid):
    """Map a SLURM job id to 'running', 'success', or 'failed'."""
    # 1. squeue first: it reflects the live controller and absorbs sacct lag.
    #    A job still listed by squeue is by definition not terminal.
    q = _run(["squeue", "-j", jobid, "-h", "-o", "%T"])
    if q is not None and q.returncode == 0 and q.stdout.strip():
        return "running"

    # 2. sacct for the terminal state. -X = allocation only, -n = no header,
    #    -P = parseable ('|' delimited, one field here).
    a = _run(["sacct", "-j", jobid, "-X", "-n", "-P", "--format=State"])
    if a is None or a.returncode != 0:
        # sacct unavailable or errored — never declare failure on no evidence;
        # let Snakemake re-poll.
        return "running"

    states = [line.strip() for line in a.stdout.splitlines() if line.strip()]
    if not states:
        # Not in the accounting DB yet (submit lag) — treat as still pending.
        return "running"

    # A state can carry a qualifier, e.g. "CANCELLED by 12345" — take the
    # leading token and normalise case.
    state = states[0].split()[0].upper()

    if state == "COMPLETED":
        return "success"
    if state in RUNNING_STATES:
        return "running"
    # Terminal failures: TIMEOUT, CANCELLED, FAILED, NODE_FAIL, OUT_OF_MEMORY,
    # BOOT_FAIL, DEADLINE, PREEMPTED, REVOKED, ...
    return "failed"


def main():
    if len(sys.argv) < 2 or not sys.argv[1].strip():
        # No job id given — safest is 'running' so Snakemake re-polls rather than
        # tearing down a job we simply could not identify.
        print("running")
        return
    # Snakemake passes what the submit command printed. sbatch --parsable emits
    # "<jobid>" (or "<jobid>;<cluster>" under multi-cluster); strip a ';cluster'
    # suffix and any '.batch'/'.extern' step suffix to leave the allocation id.
    jobid = sys.argv[1].strip().split(";")[0].split(".")[0]
    print(classify(jobid))


if __name__ == "__main__":
    main()
