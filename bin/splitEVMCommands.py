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

    results_dir = Path(os.environ.get("SYLVAN_RESULTS_DIR", "results"))
    if not results_dir.is_absolute():
        results_dir = (Path.cwd() / results_dir).resolve()

    evm_dir = Path(os.environ.get("SYLVAN_EVM_DIR", results_dir / "EVM"))
    if not evm_dir.is_absolute():
        evm_dir = (Path.cwd() / evm_dir).resolve()

    input_file = evm_dir / "commands.list"
    lines = input_file.read_text().splitlines(keepends=True)

    fixed_lines = []
    evm_path = f"{evm_dir}/"
    for line in lines:
        # Replace bare EVM/ with the resolved EVM directory; ignore occurrences already inside other paths
        fixed_lines.append(re.sub(r"(?<![\\w/])EVM/", evm_path, line))

    random.shuffle(fixed_lines)

    total_lines = len(fixed_lines)
    lines_per_file = total_lines // num_files
    remainder_lines = total_lines % num_files

    for file_num in range(num_files):
        output_file = evm_dir / f"commands.{file_num}.list"
        lines_to_write = lines_per_file + (1 if file_num < remainder_lines else 0)
        with output_file.open("w") as outfile:
            for _ in range(lines_to_write):
                outfile.write(fixed_lines.pop())


if __name__ == "__main__":
    main()
