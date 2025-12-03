#!/usr/bin/env python
"""
Generate a cluster_annotate.yml from a config_annotate.yml.

Usage:
  python bin/generate_cluster_from_config.py \
    --config config/config_annotate.yml \
    --out config/cluster_annotate.yml \
    --account cpu-s1-pgl-0 --partition cpu-s1-pgl-0

Only top-level sections that contain ncpus/threads/memory/time are copied.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict

import yaml


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Generate cluster_annotate.yml from config_annotate.yml")
    p.add_argument("--config", required=True, help="Path to config_annotate.yml")
    p.add_argument("--out", required=True, help="Output path for cluster_annotate.yml")
    p.add_argument("--account", default="placeholder", help="SLURM account (default: placeholder)")
    p.add_argument("--partition", default="placeholder", help="SLURM partition (default: placeholder)")
    p.add_argument("--time", default="14-00:00:00", help="Walltime (default: 14-00:00:00)")
    p.add_argument("--default-memory", default="4g", help="Default memory for __default__ (default: 4g)")
    p.add_argument("--default-ncpus", type=int, default=1, help="Default CPUs for __default__ (default: 1)")
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


def main() -> None:
    args = build_arg_parser().parse_args()
    cfg_path = Path(args.config)
    out_path = Path(args.out)

    cfg = load_config(cfg_path)

    cluster: Dict[str, Any] = {
        "__default__": {
            "account": args.account,
            "partition": args.partition,
            "memory": args.default_memory,
            "ncpus": args.default_ncpus,
            "nodes": 1,
            "time": args.time,
            "name": "{rule}.{wildcards}",
            "output": args.output_pattern,
            "error": args.error_pattern,
        }
    }

    for name, section in cfg.items():
        if isinstance(section, dict) and should_copy_section(section):
            cluster[name] = {}
            for key in ("ncpus", "threads", "memory", "time"):
                if key in section:
                    cluster[name][key] = section[key]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        yaml.safe_dump(cluster, fh, sort_keys=False, default_flow_style=False)

    print(f"Wrote cluster config to {out_path}")


if __name__ == "__main__":
    main()
