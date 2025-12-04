#!/usr/bin/env python3
"""Leave-one-feature-out analysis for the semi-supervised filter.

This script reuses the `Filter.semiSupRandomForest` training loop and reruns it
multiple times while removing one feature at a time. The difference in the
final out-of-bag (OOB) error provides an intuitive importance score: if dropping
an evidence track increases the OOB error, that feature contributes useful
signal to the classifier.

Example:
    python bin/filter_feature_importance.py FILTER/data.tsv results/busco/run.tsv \
        --output-table FILTER/feature_importance.tsv
"""
from __future__ import annotations

import argparse
import json
import math
import os
from typing import Dict, List

import pandas as pd

from Filter import semiSupRandomForest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Leave-one-feature-out analysis for the Sylvan filter"
    )
    parser.add_argument(
        "data",
        help="Path to the TSV created by Filter.filter_genes (e.g. FILTER/data.tsv)",
    )
    parser.add_argument(
        "busco",
        help=(
            "Path to the BUSCO table used for monitoring (same input passed to "
            "Filter.py)."
        ),
    )
    parser.add_argument(
        "--features",
        nargs="*",
        default=None,
        help=(
            "Explicit list of feature columns to evaluate. The default uses all "
            "columns except transcript_id/label and anything listed via --ignore."
        ),
    )
    parser.add_argument(
        "--ignore",
        nargs="*",
        default=[],
        help="Columns in the data TSV that should never be used as model features.",
    )
    parser.add_argument(
        "--trees",
        type=int,
        default=100,
        help="Number of trees per random forest run (default: 100)",
    )
    parser.add_argument(
        "--predictors",
        type=int,
        default=6,
        help="max_features hyperparameter for RandomForestClassifier (default: 6)",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=5,
        help=(
            "Maximum number of recycling iterations used by the semi-supervised "
            "training loop (default: 5)"
        ),
    )
    parser.add_argument(
        "--recycle",
        type=float,
        default=0.95,
        help=(
            "Prediction probability threshold required to recycle unlabeled "
            "examples (default: 0.95)"
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=123,
        help="Random seed passed to RandomForestClassifier (default: 123)",
    )
    parser.add_argument(
        "--output-table",
        default=None,
        help=(
            "Output TSV summarizing the baseline run and each leave-one-feature-out "
            "experiment (default: <data_dir>/feature_importance.tsv)"
        ),
    )
    parser.add_argument(
        "--output-json",
        default=None,
        help=(
            "Optional JSON file capturing the same summary (default: "
            "<data_dir>/feature_importance.json)"
        ),
    )
    return parser.parse_args()


def resolve_feature_list(df: pd.DataFrame, include: List[str] | None, ignore: List[str]) -> List[str]:
    """Return the ordered feature list used for training/ablation."""
    metadata_cols = {"transcript_id", "label"}
    metadata_cols.update(ignore or [])
    default_features = [c for c in df.columns if c not in metadata_cols]

    if include:
        missing = sorted(set(include) - set(default_features))
        if missing:
            raise ValueError(
                f"Requested feature(s) not found in data columns: {', '.join(missing)}"
            )
        return include

    return default_features


def summarize_process(process: Dict[str, List[float]]) -> Dict[str, float]:
    """Extract final iteration statistics from the training process log."""
    def last(seq: List[float]) -> float:
        if not seq:
            return float("nan")
        return seq[-1]

    return {
        "iterations": len(process.get("kept", [])),
        "final_kept": last(process.get("kept", [])),
        "final_discarded": last(process.get("discarded", [])),
        "final_kept_buscos": last(process.get("kept_buscos", [])),
        "final_discarded_buscos": last(process.get("discarded_buscos", [])),
        "final_oob_error": last(process.get("OOB", [])),
    }


def run_filter(
    data: pd.DataFrame,
    features: List[str],
    busco_path: str,
    args: argparse.Namespace,
) -> Dict[str, float]:
    subset_cols = ["transcript_id", "label"] + features
    subset = data.loc[:, subset_cols].copy()
    _, process = semiSupRandomForest(
        subset,
        args.predictors,
        busco_path,
        args.trees,
        seed=args.seed,
        recycle_prob=args.recycle,
        maxiter=args.max_iter,
    )
    return summarize_process(process)


def format_delta(value: float) -> str:
    if value is None or math.isnan(value):
        return "nan"
    return f"{value:+.4f}"


def main() -> None:
    args = parse_args()
    data = pd.read_csv(args.data, sep="\t")

    feature_list = resolve_feature_list(data, args.features, args.ignore)
    if not feature_list:
        raise ValueError("No usable features detected in data TSV.")

    out_dir = os.path.dirname(os.path.abspath(args.data))
    table_path = args.output_table or os.path.join(out_dir, "feature_importance.tsv")
    json_path = args.output_json or os.path.join(out_dir, "feature_importance.json")

    print(f"Running baseline model with {len(feature_list)} features...")
    baseline = run_filter(data, feature_list, args.busco, args)
    baseline_row = {
        "feature_removed": "(none)",
        "num_features": len(feature_list),
        "oob_delta": 0.0,
        **baseline,
    }

    results = [baseline_row]
    for feature in feature_list:
        reduced = [f for f in feature_list if f != feature]
        if not reduced:
            continue
        print(f"Dropping '{feature}' ({len(reduced)} features remaining)...")
        summary = run_filter(data, reduced, args.busco, args)
        summary_row = {
            "feature_removed": feature,
            "num_features": len(reduced),
            "oob_delta": summary["final_oob_error"] - baseline["final_oob_error"],
            **summary,
        }
        results.append(summary_row)
        delta_str = format_delta(summary_row["oob_delta"])
        print(
            f"  -> final OOB error: {summary_row['final_oob_error']:.4f} "
            f"(delta {delta_str})"
        )

    df = pd.DataFrame(results)
    df.to_csv(table_path, sep="\t", index=False)
    print(f"\nSummary written to {table_path}")

    if json_path:
        json_ready = []
        for row in results:
            json_ready.append(
                {
                    k: (None if isinstance(v, float) and math.isnan(v) else v)
                    for k, v in row.items()
                }
            )
        with open(json_path, "w", encoding="utf-8") as fh:
            json.dump({"runs": json_ready}, fh, indent=2)
        print(f"JSON summary written to {json_path}")


if __name__ == "__main__":
    main()
