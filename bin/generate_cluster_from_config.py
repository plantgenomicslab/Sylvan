#!/usr/bin/env python
"""
Generate a standalone cluster config YAML from a config_annotate.yml.

Note: config_annotate.yml can be used directly as --cluster-config for Snakemake
(extra top-level keys like Input and Internal are ignored). This script is only
needed if you want a separate, minimal cluster-only file.

Usage:
  python bin/generate_cluster_from_config.py \
    --config config/config_annotate.yml \
    --out cluster_annotate.yml \
    --account cpu-s1-pgl-0 --partition cpu-s1-pgl-0

Only top-level sections that contain ncpus/threads/memory/time are copied.
"""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path
from typing import Any, Dict

import yaml

FALLBACK_TIME = "9-00:00:00"


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generate a standalone cluster YAML from config_annotate.yml"
    )
    p.add_argument("--config", required=True, help="Path to config_annotate.yml")
    p.add_argument("--out", required=True, help="Output path for cluster YAML")
    p.add_argument(
        "--account", default="placeholder", help="SLURM account (default: placeholder)"
    )
    p.add_argument(
        "--partition",
        default="placeholder",
        help="SLURM partition (default: placeholder)",
    )
    p.add_argument(
        "--time", default="auto",
        help="Walltime (default: auto = partition max - 1 day, fallback %(default)s)"
    )
    p.add_argument(
        "--default-memory",
        default="4g",
        help="Default memory for __default__ (default: 4g)",
    )
    p.add_argument(
        "--default-ncpus",
        type=int,
        default=1,
        help="Default CPUs for __default__ (default: 1)",
    )
    p.add_argument(
        "--output-pattern",
        default="results/logs/{rule}_{wildcards}.out",
        help="Output log pattern for __default__",
    )
    p.add_argument(
        "--error-pattern",
        default="results/logs/{rule}_{wildcards}.err",
        help="Error log pattern for __default__",
    )
    return p


def load_config(path: Path) -> Dict[str, Any]:
    with path.open() as fh:
        return yaml.safe_load(fh) or {}


def should_copy_section(section: Dict[str, Any]) -> bool:
    return any(key in section for key in ("ncpus", "threads", "memory", "time"))


def _parse_slurm_time_to_seconds(time_str):
    """Parse a SLURM time string (D-HH:MM:SS, HH:MM:SS, etc.) to total seconds."""
    days = 0
    if "-" in time_str:
        day_part, time_str = time_str.split("-", 1)
        days = int(day_part)
    parts = time_str.split(":")
    parts = [int(p) for p in parts]
    if len(parts) == 3:
        hours, minutes, seconds = parts
    elif len(parts) == 2:
        hours, minutes = parts
        seconds = 0
    else:
        hours = parts[0]
        minutes = seconds = 0
    return days * 86400 + hours * 3600 + minutes * 60 + seconds


def _get_partition_max_time(partition):
    """Query sinfo for the partition's max time limit. Returns D-HH:MM:SS or None."""
    if not partition or partition.lower() in ("", "placeholder"):
        return None
    try:
        result = subprocess.run(
            ["sinfo", "-p", partition, "-h", "-o", "%l"],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode != 0 or not result.stdout.strip():
            return None
        max_time = result.stdout.strip().splitlines()[0].strip()
        if max_time.lower() == "infinite":
            return None
        return max_time
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None


def _auto_time(partition):
    """Return (max_time - 1 day) for the partition. Falls back to FALLBACK_TIME on failure."""
    max_time_str = _get_partition_max_time(partition)
    if not max_time_str:
        print(f"WARNING: Could not read max time from sinfo, using fallback: {FALLBACK_TIME}")
        return FALLBACK_TIME
    total_seconds = _parse_slurm_time_to_seconds(max_time_str)
    adjusted = total_seconds - 86400  # subtract 1 day
    if adjusted <= 0:
        print(f"WARNING: Partition max time ({max_time_str}) is <= 1 day, using fallback: {FALLBACK_TIME}")
        return FALLBACK_TIME
    return _normalize_time(adjusted)


def _normalize_time(val):
    """Convert a time value (possibly parsed as sexagesimal int) to D-HH:MM:SS string.

    YAML 1.1 parses unquoted colon-separated values like 72:00:00 as sexagesimal
    integers (72*3600 + 0*60 + 0 = 259200).  This converts them back to a valid
    SLURM walltime string.
    """
    if isinstance(val, str):
        return val
    total_seconds = int(val)
    days, remainder = divmod(total_seconds, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{days}-{hours:02d}:{minutes:02d}:{seconds:02d}"


def main() -> None:
    args = build_arg_parser().parse_args()
    cfg_path = Path(args.config)
    out_path = Path(args.out)

    cfg = load_config(cfg_path)

    # Resolve auto time: query sinfo for partition max - 1 day
    time_value = args.time
    if time_value.lower() == "auto":
        time_value = _auto_time(args.partition)
        print(f"Auto time: {time_value} (partition: {args.partition})")

    cluster: Dict[str, Any] = {
        "__default__": {
            "account": args.account,
            "partition": args.partition,
            "memory": args.default_memory,
            "ncpus": args.default_ncpus,
            "nodes": 1,
            "time": time_value,
            "name": "{rule}.{wildcards}",
            "output": args.output_pattern,
            "error": args.error_pattern,
            "extra_args": "",
        }
    }

    for name, section in cfg.items():
        if name == "__default__":
            continue  # Already built from CLI args above
        if isinstance(section, dict) and should_copy_section(section):
            cluster[name] = {}
            for key in ("ncpus", "threads", "memory", "time", "extra_args"):
                if key in section:
                    val = section[key]
                    # YAML 1.1 parses unquoted colon-separated values (e.g. 72:00:00)
                    # as sexagesimal integers.  Coerce time back to a string so sbatch
                    # receives a valid walltime instead of a huge minute count.
                    if key == "time":
                        val = _normalize_time(val)
                    cluster[name][key] = val

    out_path.parent.mkdir(parents=True, exist_ok=True)
    raw_yaml = yaml.safe_dump(cluster, sort_keys=False, default_flow_style=False)
    # Ensure time values are always quoted to prevent YAML 1.1 sexagesimal parsing
    raw_yaml = re.sub(
        r'^(\s*time:\s*)(\S.*)$',
        lambda m: m.group(1) + (m.group(2) if m.group(2).startswith(("'", '"')) else f'"{m.group(2)}"'),
        raw_yaml,
        flags=re.MULTILINE,
    )
    with out_path.open("w") as fh:
        fh.write(raw_yaml)

    print(f"Wrote cluster config to {out_path}")


if __name__ == "__main__":
    main()
