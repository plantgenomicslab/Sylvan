#!/usr/bin/env python
import sys
import math
import random
import os
import re

def split_file(input_file, num_files):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # Fix paths in commands: replace standalone EVM/ with results/EVM/
    # This handles cases where partitions_list.out was generated with old paths
    fixed_lines = []
    for line in lines:
        # Replace --exec_dir EVM/ with --exec_dir results/EVM/
        line = re.sub(r'--exec_dir\s+EVM/', '--exec_dir results/EVM/', line)
        # Replace > EVM/ with > results/EVM/ (for output redirection)
        line = re.sub(r'>\s*EVM/', '> results/EVM/', line)
        # Replace 2> EVM/ with 2> results/EVM/ (for stderr redirection)
        line = re.sub(r'2>\s*EVM/', '2> results/EVM/', line)
        fixed_lines.append(line)

    random.shuffle(fixed_lines)  # Shuffle the lines

    total_lines = len(fixed_lines)
    lines_per_file = total_lines // num_files
    remainder_lines = total_lines % num_files

    for file_num in range(num_files):
        output_file = f"results/EVM/commands.{file_num}.list"
        with open(output_file, 'w') as outfile:
            lines_to_write = lines_per_file
            if file_num < remainder_lines:
                lines_to_write += 1
            for _ in range(lines_to_write):
                outfile.write(fixed_lines.pop())

# Example usage
input_file = "results/EVM/commands.list"  # Specify your input file here
num_files = int(sys.argv[1])  # Specify the desired number of output files
split_file(input_file, num_files)

