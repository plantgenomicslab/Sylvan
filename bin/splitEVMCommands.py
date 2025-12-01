#!/usr/bin/env python
import os
import random
import re
import sys
from pathlib import Path


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: splitEVMCommands.py <num_files>")

    num_files = int(sys.argv[1])

    # Default to repo-root/results/EVM; allow env overrides that can be absolute or relative to CWD
    root_dir = Path(__file__).resolve().parent.parent
    default_results = root_dir / "results"

    results_env = os.environ.get("SYLVAN_RESULTS_DIR")
    if results_env:
        results_dir = Path(results_env)
        if not results_dir.is_absolute():
            results_dir = (Path.cwd() / results_env).resolve()
    else:
        results_dir = default_results.resolve()

    evm_env = os.environ.get("SYLVAN_EVM_DIR")
    if evm_env:
        evm_dir = Path(evm_env)
        if not evm_dir.is_absolute():
            evm_dir = (Path.cwd() / evm_env).resolve()
    else:
        evm_dir = (results_dir / "EVM").resolve()

    evm_dir.mkdir(parents=True, exist_ok=True)

    input_file = evm_dir / "commands.list"
    lines = input_file.read_text().splitlines(keepends=True)

    fixed_lines = []
    for line in lines:
        # The evidence_modeler.pl script changes its working directory to the value of --exec_dir.
        # This means that any relative paths for input files will be broken.
        # This loop finds common input file arguments and converts their paths to be absolute.
        # It looks for a flag (e.g., -G), whitespace, and then a path that does *not* start with '/'.
        # It assumes file paths do not contain spaces.
        updated_line = line
        for flag in ["-G", "-g", "-e", "-p"]:
            # Pattern: flag, one or more spaces, a path that is not absolute (doesn't start with /),
            # and is followed by a space. We'll capture the flag and path.
            pattern = re.compile(f"({flag}\\s+)([^/\\s][^\\s]*)")

            # Use a function for replacement to handle multiple matches in a line if they were to occur
            def repl(match):
                file_path = match.group(2)
                # Prepend the absolute path to the EVM directory
                abs_path = evm_dir / file_path
                return f"{match.group(1)}{abs_path}"

            updated_line = pattern.sub(repl, updated_line)

        fixed_lines.append(updated_line)

    random.shuffle(fixed_lines)

    total_lines = len(fixed_lines)
    lines_per_file = total_lines // num_files
    remainder_lines = total_lines % num_files

    for file_num in range(num_files):
        output_file = evm_dir / f"commands.{file_num}.list"
        lines_to_write = lines_per_file + (1 if file_num < remainder_lines else 0)
        with output_file.open("w") as outfile:
            for _ in range(lines_to_write):
                if fixed_lines:
                    outfile.write(fixed_lines.pop())


if __name__ == "__main__":
    main()
