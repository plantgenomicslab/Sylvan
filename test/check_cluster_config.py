#!/usr/bin/env python3
"""Compare a Snakefile's rules against its cluster-config resource keys.

Snakemake matches `--cluster-config` blocks by rule name. A key that does not
name a rule is inert, and a rule with no key silently takes `__default__` --
which is how issue #14 got 16-thread jobs onto a 4-cpu allocation without
anything failing at submission time.

Three checks, in the order they cost real time:

1. threads > ncpus   -- oversubscription. This is the one that keeps finding
                        live defects: the right-sizing pass lowered
                        calculateRSEMExpression_PAIRED to 12 cpus but left the
                        shared thread key at 16.
2. inert key         -- a config block naming no rule in this workflow.
3. rule on default   -- reported, not failed: some rules take __default__ on
                        purpose.

Cross-workflow keys are called out separately. config_filter.yml carried
`extractTranscriptsForFilter` and `filterBlastP`, which are annotate rules --
inert here, and disagreeing with the annotate values that do apply.
"""

import argparse
import re
import sys
from pathlib import Path

import yaml

# Top-level blocks that configure the pipeline rather than a rule.
RESERVED = {
    "Input", "Cutoff", "Threads", "Benchmark",
    "Singularity", "Internal", "__default__",
}

RULE_RE = re.compile(r"(?m)^\s*(?:rule|checkpoint)\s+([A-Za-z_][A-Za-z0-9_]*)\s*:")
THREADS_KEY_RE = re.compile(r"""Threads["']\]\[["']([A-Za-z_][A-Za-z0-9_]*)["']\]""")


def rule_names(snakefile):
    return set(RULE_RE.findall(Path(snakefile).read_text(encoding="utf-8")))


def rule_threads(snakefile, threads_section):
    """rule -> declared thread count, for rules that declare one resolvably.

    Handles `threads: 8` and `threads: config_dict["Threads"]["key"]`. A
    `.get(...)` fallback or any other expression is skipped rather than guessed
    at -- an unresolved expression is not evidence of a defect.
    """
    text = Path(snakefile).read_text(encoding="utf-8")
    out = {}
    starts = [(m.start(), m.group(1)) for m in RULE_RE.finditer(text)]
    for i, (pos, name) in enumerate(starts):
        end = starts[i + 1][0] if i + 1 < len(starts) else len(text)
        block = text[pos:end]
        m = re.search(r"(?m)^\s*threads:\s*(.+)$", block)
        if not m:
            continue
        expr = m.group(1).strip()
        if expr.isdigit():
            out[name] = int(expr)
            continue
        km = THREADS_KEY_RE.search(expr)
        if km:
            val = threads_section.get(km.group(1))
            if isinstance(val, int):
                out[name] = val
            elif val is None:
                out[name] = f"MISSING:{km.group(1)}"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snakefile", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--threads-from",
                    help="config carrying the Threads section, when the cluster "
                         "resources live in a separate file (toydata splits them). "
                         "Without a Threads section the thread check is skipped "
                         "rather than reported as missing.")
    ap.add_argument("--also-rules-from", action="append", default=[],
                    help="Snakefile whose rules make a key cross-workflow "
                         "rather than unknown (reported, still inert).")
    args = ap.parse_args()

    config = yaml.safe_load(Path(args.config).read_text(encoding="utf-8")) or {}
    keys = {k for k, v in config.items()
            if k not in RESERVED and isinstance(v, dict)}
    rules = rule_names(args.snakefile)
    other = set()
    for extra in args.also_rules_from:
        other |= rule_names(extra)

    default_cpus = (config.get("__default__") or {}).get("ncpus")
    threads_src = config
    if args.threads_from:
        threads_src = yaml.safe_load(
            Path(args.threads_from).read_text(encoding="utf-8")) or {}
    threads_section = threads_src.get("Threads") or {}
    threads = rule_threads(args.snakefile, threads_section) if threads_section else {}
    if not threads_section:
        print("note: no Threads section in "
              f"{args.threads_from or args.config}; skipping the threads check")

    errors = []
    warnings = []

    # 1. oversubscription
    for rule in sorted(rules):
        declared = threads.get(rule)
        if isinstance(declared, str) and declared.startswith("MISSING:"):
            errors.append(f"{rule}: threads reads Threads.{declared[8:]}, "
                          f"which is not defined")
            continue
        if not isinstance(declared, int):
            continue
        ncpus = (config.get(rule) or {}).get("ncpus", default_cpus)
        if isinstance(ncpus, int) and declared > ncpus:
            errors.append(f"{rule}: threads {declared} > ncpus {ncpus} "
                          f"({'__default__' if rule not in keys else 'own key'})")

    # 2. inert keys
    for key in sorted(keys - rules):
        where = " (a rule of another workflow)" if key in other else ""
        errors.append(f"{key}: config block names no rule in {args.snakefile}{where}")

    # 3. rules on __default__
    for rule in sorted(rules - keys):
        if rule == "all":
            continue
        warnings.append(f"{rule}: no key, takes __default__")

    for line in warnings:
        print(f"note: {line}")
    for line in errors:
        print(f"error: {line}")

    print(f"\n{len(rules)} rules, {len(keys)} keys, "
          f"{len(errors)} error(s), {len(warnings)} note(s)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
