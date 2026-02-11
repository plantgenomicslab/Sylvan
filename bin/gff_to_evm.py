#!/usr/bin/env python 
# Convert input files to EVM format
# Usage: gff_to_evm.py [file] [identifier] > [outfile]
#
# Note: identifier matches evm config (e.g. nucleotide_to_protein_match)

import re
import sys
import argparse

ap = argparse.ArgumentParser(
	prog='gff_to_evm.py',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	description=('Reformat evidence files for evm.')
)
ap.add_argument('file', type=str, help="Input file to be reformatted")
ap.add_argument('--source', type=str, default="", help="Evidence source (matches entry in EVM weights file)")
ap.add_argument('--type', type=str, default="", help="Evidence type (goes in gff column 3)")
ap.add_argument('--check-coords', action='store_true', help="If set will ensure start coordinate is less than end coordinate")
ap.add_argument('--miniprot', action='store_true', help="Attribute handling for miniprot")
ap.add_argument('--genewise', action='store_true', help="Attribute handling for genewise")
ap.add_argument('--feature', type=str, default="CDS", help="Feature type to extract. (GFF column 3 )")
args = ap.parse_args()

with open(args.file, "r") as f:
	for line in f:
		if not line.startswith("#") and re.search(f"\t{args.feature}\t", line):
			if args.genewise:
				line = re.sub("Name=", "Target=", line)
			if args.miniprot:
				line = re.sub("Parent=", "ID=", line)
		
			fields = line.strip().split("\t")
			
			if args.source != "":
				fields[1] = args.source
			if args.type != "":
				fields[2] = args.type
		
			if args.check_coords:
				if int(fields[3]) > int(fields[4]):
					start = fields[4]
					end = fields[3]
					fields[3] = start
					fields[4] = end
					#fields[3], fields[4] = fields[4], fields[3]

			print("\t".join(fields[:9]))
