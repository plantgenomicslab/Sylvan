#!/usr/bin/env python
import sys
import math
import random
import os

def split_file(input_file, num_files):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()
    random.shuffle(lines)  # Shuffle the lines

    total_lines = len(lines)
    lines_per_file = total_lines // num_files
    remainder_lines = total_lines % num_files

    for file_num in range(num_files):
        output_file = f"results/EVM/commands.{file_num}.list"
        with open(output_file, 'w') as outfile:
            lines_to_write = lines_per_file
            if file_num < remainder_lines:
                lines_to_write += 1
            for _ in range(lines_to_write):
                outfile.write(lines.pop())

# Example usage
input_file = "results/EVM/commands.list"  # Specify your input file here
num_files = int(sys.argv[1])  # Specify the desired number of output files
split_file(input_file, num_files)

