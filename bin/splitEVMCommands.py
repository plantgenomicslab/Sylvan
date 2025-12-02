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

    # Get the absolute path to the EVM directory from environment variables
    # These are set in the Snakefile just before this script is called.
    evm_env = os.environ.get("SYLVAN_EVM_DIR")
    if not evm_env:
        sys.exit("Error: SYLVAN_EVM_DIR environment variable not set.")
    
    # Resolve to be absolutely sure it's an absolute path
    evm_dir = Path(evm_env).resolve()

    input_file = evm_dir / "commands.list"
    if not input_file.is_file():
        # If the file doesn't exist, there are no commands to run.
        # Create empty command files to satisfy Snakemake dependencies.
        for file_num in range(num_files):
            output_file = evm_dir / f"commands.{file_num}.list"
            output_file.touch()
        return
        
    lines = input_file.read_text().splitlines(keepends=True)

    fixed_lines = []
    for line in lines:
        updated_line = line
        # These are the flags that take a single file path argument that needs to be made absolute.
        for flag in ["-G", "-g", "-e", "-p"]:
            # The pattern looks for:
            # - the flag (e.g., "-G") followed by one or more whitespace characters.
            # - It captures the flag and whitespace in group 1.
            # - It then captures a sequence of characters that is the file path in group 2.
            #   - The path is assumed not to start with a "/" (i.e., it's a relative path).
            #   - The path consists of any characters that are not whitespace.
            pattern = re.compile(f"({flag}\\s+)([^/\\s][^\\s]*)")
            
            # The replacement function constructs the absolute path.
            def repl(match):
                file_path = match.group(2)
                # All relative paths in commands.list are relative to the EVM directory.
                abs_path = evm_dir / file_path
                return f"{match.group(1)}{abs_path}"

            updated_line = pattern.sub(repl, updated_line)
        
        fixed_lines.append(updated_line)

    # Shuffle and distribute the fixed commands
    random.shuffle(fixed_lines)

    total_lines = len(fixed_lines)
    if total_lines == 0:
        # Create empty files if after processing there are no commands
        for file_num in range(num_files):
             output_file = evm_dir / f"commands.{file_num}.list"
             output_file.touch()
        return

    # Integer division gives the base number of lines per file
    lines_per_file = total_lines // num_files
    # The remainder is distributed one by one to the first few files
    remainder_lines = total_lines % num_files

    current_line_idx = 0
    for file_num in range(num_files):
        output_file = evm_dir / f"commands.{file_num}.list"
        num_to_write = lines_per_file + (1 if file_num < remainder_lines else 0)
        
        with output_file.open("w") as outfile:
            for _ in range(num_to_write):
                if current_line_idx < total_lines:
                    outfile.write(fixed_lines[current_line_idx])
                    current_line_idx += 1

if __name__ == "__main__":
    main()
