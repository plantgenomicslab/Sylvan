#!/usr/bin/env python
"""
https://github.com/plantgenomicslab/geta/blob/master/bin/remove_genes_in_repeats
"""
import argparse
import sys
import re
import pandas as pd

parser = argparse.ArgumentParser(description='Remove gene models that overlap with repeats in a genome GFF3 file')
parser.add_argument('repeat_gff3', type=str, help='GFF3 file of repeats')
parser.add_argument('genome_gff3', type=str, help='GFF3 file of genome')
parser.add_argument('--ratio', type=float, default=0.8, help='The ratio of CDS_overlap_region_base_number / CDS_total_base_number')
parser.add_argument('--filtered_gene_models', type=str, default=None, help='Output the filtered gene Models in GFF3 format')
parser.add_argument('--ignore_Simple_repeat', action='store_true', help='Ignore Simple_repeat and Satellite repeat types')
parser.add_argument('--ignore_Unknown', action='store_true', help='Ignore Unknown repeat types')
args = parser.parse_args()

repeat_df = pd.read_csv(args.repeat_gff3, sep='\t', comment='#', header=None)
genome_df = pd.read_csv(args.genome_gff3, sep='\t', comment='#', header=None)


def cal_overlap_length(repeat, i: str) -> int:
    i = i.split("\t")
    length = 0
    index1 = int(int(i[1]) / 1000)
    index2 = int(int(i[2]) / 1000)
    repeat_info = {}
    for j in range(index1, index2+1):
        if i[0] in repeat:
            if j in repeat[i[0]]:
                for k in repeat[i[0]][j]:
                    repeat_info[k] = 1
    repeat_info_keys = list(repeat_info.keys())
    repeat_info_keys.sort()
    
    start = int(i[1])
    for rr in repeat_info_keys:
        r = rr.split("\t")
        if int(r[0]) <= int(i[2]) and int(r[1]) >= start:
            if start <= int(r[0]):
                if int(i[2]) < int(r[1]):
                    length += (int(i[2]) - int(r[0]) + 1)
                    break
                else:
                    length += (int(r[1]) - int(r[0]) + 1)
                    start = int(r[1]) + 1
            else:
                if int(i[2]) < int(r[1]):
                    length += (int(i[2]) - start + 1)
                    break
                else:
                    length += (int(r[1]) - start + 1)
                    start = int(r[1]) + 1
    
    return length
    
   # logic to remove gene models that overlap with repeats
    # ...

# Get repeat file name from command line argument
repeat_file = sys.argv[1]
repeat = {}

# Exit the script if the file cannot be opened
try:
    with open(repeat_file, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip comment lines and empty lines
            if line.startswith("#") or not line:
                continue

            # Skip Simple_repeat and Satellite repeat types if ignore_Simple_repeat argument passed
            if args.ignore_Simple_repeat:
                if "Name=" in line:
                    repeat_type = line.split("Name=")[1].split(";")[0]
                    if repeat_type in ["Simple_repeat", "Satellite"]:
                        continue

            # Skip Unknown repeat types if ignore_Unknown argument passed
            if args.ignore_Unknown:
                if "Name=" in line:
                    repeat_type = line.split("Name=")[1].split(";")[0]
                    if repeat_type == "Unknown":
                        continue

            fields = line.split("\t")
            # print(fields[5])
            # print(fields[3])
            # print(fields[4])
            # print(fields[0])
            index1 = int(int(fields[3]) / 1000)
            index2 = int(int(fields[4]) / 1000)
            for i in range(index1, index2 + 1):
                if fields[0] not in repeat:
                    repeat[fields[0]] = {}
                if i not in repeat[fields[0]]:
                    repeat[fields[0]][i] = {}
                repeat[fields[0]][i][f"{fields[3]}\t{fields[4]}"] = 1
except OSError:
    print("Cannot open file:", repeat_file)
    sys.exit()

# Get genome file name from command line argument
genome_file = sys.argv[2]

# Exit the script if the file cannot be opened
try:
    with open(genome_file, 'r') as f:
        gene = {}
        gene_id = None
        gene_id_list = []
        for line in f:
            match = re.search(r"\tgene\t.*ID=([^\s;]+)", line)
            if match:
                gene_id = match.group(1)
                gene_id_list.append(gene_id)
            if gene_id:
                gene[gene_id] = gene.get(gene_id, "") + line
except OSError:
    print("Cannot open file:", genome_file)
    sys.exit()

for gene_id in gene_id_list:
    lines = gene[gene_id].split("\n")
    cds = [line for line in lines if "\tCDS\t" in line]
    if len(cds) == 0:
        continue

    total_length = 0
    overlap_length = 0
    for c in cds:
        fields = c.split("\t")
        #print({fields[0]},"\t",{fields[3]},"\t",{fields[4]},"\t",{fields[5]})
        overlap_length += cal_overlap_length(repeat, f"{fields[0]}\t{fields[3]}\t{fields[4]}")
        #print(overlap_length, file=sys.stderr)
        total_length += (int(fields[4]) - int(fields[3]) + 1)

    coverage = overlap_length / total_length
    if coverage >= args.ratio:
        print(f"{gene_id}\t{overlap_length}\t{total_length}", file=sys.stderr)
        if args.filtered_gene_models:
            with open(args.filtered_gene_models, 'a') as out:
                out.write(gene[gene_id])
    else:
        print(gene[gene_id])
