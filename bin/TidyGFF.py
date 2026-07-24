import argparse
import re
import math
import os
import subprocess
import pandas as pd

def getParent(attr: str) -> str:
	# Anchor to start/';' and accept any non-';'/space char (issue #20.9): the old
	# [a-zA-Z\d.\-_:] class truncated IDs containing '|' etc. and mismatched
	# replaceParent's class.
	m = re.search(r"(?:^|;)Parent=([^;\s]+)", attr)
	if not m:
		raise ValueError(f"No Parent= found in attribute string: {attr}")
	return m.group(1)

def getID(attr:str) -> str:
	m = re.search(r"(?:^|;)ID=([^;\s]+)", attr)
	if not m:
		raise ValueError(f"No ID= found in attribute string: {attr}")
	return m.group(1)

def findChroms(chrom: str, chrom_regex = False) -> str:
	pre = ""
	if chrom_regex:
		if re.search(chrom_regex, chrom):
			pre = re.search(chrom_regex, chrom).group(0)
	elif re.search(r"(^Chr)|(^chr)|(^LG)|(^Ch)|(^\d)", chrom):
		pre = re.search(r"(^Chr)|(^chr)|(^LG)|(^Ch)|(^\d)", chrom).group(0)
	# No prefix matched (e.g. ptg/scaffold contigs not covered by the regex):
	# return the name unchanged instead of crashing on an undefined 'pre'.
	if not pre:
		return chrom
	return re.sub(pre, "", chrom)

def replaceParent(id: str, attr: str) -> str:
	replace = f"Parent={id};"
	# count=1: only the real Parent= is rewritten, never a 'Parent=' substring
	# inside a Note. Broadened char class matches getParent (issue #20.9).
	return re.sub(r"Parent=[^;\s]*;?", replace, attr, count=1)

def replaceID(id: str, attr: str, suf: str) -> str:
	replace = f"ID={id}{suf};"
	replace = re.sub(r"ID=[^;\s]*;?", replace, attr, count=1)
	old_m = re.search(r"ID=([^;\s]*);?", attr)
	old = old_m.group(1) if old_m else ""
	return replace, old

def printData(out, item, children):
	data = [item["data"][entry] for entry in item["data"].keys()]
	out.write("\t".join(data) + "\n")

	if children:
		children = [i for i in item.keys() if i != "data"]

		if len(children) > 0:
			for i in children:
				child = item[i]
				printData(out, child, True)

def upgradeChromosome(chrom: str, just: int, chrom_regex = False, contig_regex = False) -> str:
	if chrom_regex:
		if re.search(chrom_regex, chrom):
			m = re.search(chrom_regex, chrom)
			chrom_num = chrom.replace(m[0], "")
			chrom_phase = re.search("[A-Z]", chrom_num)
			if chrom_phase:
				chrom_num = chrom_num.replace(chrom_phase[0], "")
				chrom_phase = chrom_phase[0]
			else:
				chrom_phase = "G"
			new_name = m[0] + str(chrom_num).zfill(just) + chrom_phase 
			return(new_name)
	elif re.search("(^Chr)|(^chr)|(^LG)|(^Ch)", chrom):
		m = re.search("(^Chr)|(^chr)|(^LG)|(^Ch)", chrom)
		remainder = chrom.replace(m[0], "", 1)
		# Fallback for non-numeric chromosomes (organellar ChrM/ChrC, suffixes):
		# int("M") would crash (issue #20.11); keep the remainder as-is instead.
		if remainder.isdigit():
			new_name = str(int(remainder)).zfill(just) + "G"
		else:
			new_name = (remainder or chrom) + "G"
		return(new_name)
	elif re.search(r"^\d", chrom):
		new_name = chrom.zfill(just) + "G"
		return(new_name)

	if contig_regex:
		if re.search(contig_regex[0], chrom):
			contig_num = re.search(r"[\d\_\-]+", chrom)[0]
			new_name = contig_regex[1] + str(contig_num) + "G"
			return(new_name)
	
	return(chrom)

def getSuffix(feature_id: str) -> str:
	feature_id = feature_id.lower()
	suf = ["cds", "exon", "intron", "utr5p", "utr3p", "five_prime_utr", "three_prime_utr"]
	for suffix in suf:
		if re.search(suffix, feature_id):
			return(re.search(rf"{suffix}\.{{0,1}}\d*", feature_id)[0].lower())

def loadGFF(gff: str) -> dict:
	with open(gff, 'r') as infile:
		gff_dict = {}
		# Resolve a feature to its mRNA's gene via a recorded map instead of the
		# last-seen gene, so an interleaved GFF does not attach features to the
		# wrong parent (issue #20.8). None => mRNA stored at chromosome top level.
		mrna_to_gene = {}
		protein_coding = True

		for line in infile:
			if line.startswith("#") or not line.strip():
				continue

			line = line.strip().split("\t")
			# Normalise field count (issue #20.8): a literal TAB inside an
			# attribute value (seen in the wild, e.g. Spe Note=) splits column 9
			# into several fields and silently truncates it. Skip lines with too
			# few fields; rejoin any surplus into column 9 with a space.
			if len(line) < 9:
				continue
			if len(line) > 9:
				line = line[:8] + [" ".join(line[8:])]

			#### For loading mikado output ###
			if (line[2] == "superlocus"):
				continue

			if (line[2] == "ncRNA_gene"):
				protein_coding = False
			##################################

			chrom = line[0]
			feature_id = getID(line[8])
			feature_data = {"seqid":chrom,
							"source":line[1],
							"type":line[2],
							"start":line[3],
							"end":line[4],
							"score": line[5],
							"strand":line[6],
							"phase":line[7],
							"attributes":line[8]}

			if chrom not in gff_dict.keys():
				gff_dict[chrom] = {}

			if line[2] == "gene":
				gff_dict[chrom][feature_id] = {"data":feature_data}
				protein_coding = True # Support for mikado
			elif line[2] == "mRNA":
				gene = getParent(line[8])

				if gene not in gff_dict[chrom].keys():
					print(f"WARNING: mRNA with id {feature_id} missing parent {gene}")
					gff_dict[chrom][feature_id] = {"data":feature_data}
					mrna_to_gene[(chrom, feature_id)] = None
				else:
					gff_dict[chrom][gene][feature_id] = {"data":feature_data}
					mrna_to_gene[(chrom, feature_id)] = gene
			else:
				if not protein_coding: #Support for mikado
					continue
				mrna = getParent(line[8])
				if (chrom, mrna) not in mrna_to_gene:
					print(f"WARNING: feature {feature_id} references unknown parent mRNA {mrna}; skipping")
					continue
				parent_gene = mrna_to_gene[(chrom, mrna)]
				if parent_gene is None:
					parent_dict = gff_dict[chrom][mrna]
				else:
					parent_dict = gff_dict[chrom][parent_gene][mrna]
				# Handle duplicate feature IDs (e.g., CDS segments sharing one ID)
				orig_id = feature_id
				suffix = 2
				while feature_id in parent_dict:
					feature_id = f"{orig_id}_{suffix}"
					suffix += 1
				parent_dict[feature_id] = {"data":feature_data}
	return(gff_dict)

def singleExonGenes(gff_dict):
	# Count EXON features, not CDS (issue #20.1). A transcript with a single CDS
	# but additional non-coding (UTR) exons is multi-exon; counting CDS mislabels
	# it single-exon and biases the RF seed labels. Fall back to CDS count only
	# when a transcript has no exon features at all (CDS-only GFF).
	singleExons = []
	for chrom in gff_dict:
		for gene in gff_dict[chrom]:
			mrna_ids = [i for i in gff_dict[chrom][gene].keys() if i != "data"]
			for mrna in mrna_ids:
				exon = 0
				cds = 0
				feature_ids = [i for i in gff_dict[chrom][gene][mrna].keys() if i != "data"]
				for feature in feature_ids:
					ftype = gff_dict[chrom][gene][mrna][feature]["data"]["type"]
					if ftype == "exon":
						exon += 1
					elif ftype == "CDS":
						cds += 1
				count = exon if exon > 0 else cds
				if count <= 1:
					singleExons.append(mrna)
	return(pd.DataFrame(singleExons, columns=["transcript_id"]))
			
def sortGFF(gff_dict: dict, output: str, chrom_regex = False, contig_regex = False) -> None:
	outfile = open(output, 'w')

	if not chrom_regex:
		chrom_regex = r"(^Chr)|(^chr)|(^LG)|(^Ch)|(^\d)"
	if not contig_regex:
		contig_regex = r"(?! (^Chr)|(^chr)|(^LG)|(^Ch)|(^\d))"
	else: 
		contig_regex = contig_regex.split(",")[0]
	
	chromSort_df = pd.DataFrame(gff_dict.keys())
	chromSort_df[1] = chromSort_df[0].apply(lambda x: findChroms(x, chrom_regex))
	chromSort_df[2] = chromSort_df[0].apply(lambda x: re.search(chrom_regex, x)[0] if re.search(chrom_regex, x) else None)
	chromSort_df[3] = chromSort_df[0].apply(lambda x: findChroms(x, contig_regex))
	chromSort_df[3] = pd.to_numeric(chromSort_df[3], errors="coerce")
	chromSort_df[4] = chromSort_df[2].apply(lambda x: "achromosome" if x != None else "bcontig")
	chromSort_chroms = chromSort_df[chromSort_df[4] == "achromosome"]
	chromSort_contigs = chromSort_df[chromSort_df[4] == "bcontig"]
	chromSort_chroms = chromSort_chroms.sort_values(by=[2,1])
	chromSort_contigs = chromSort_contigs.sort_values(by=[3])
	chromSort_df = pd.concat([chromSort_chroms,chromSort_contigs])

	for chrom in chromSort_df[0]:
		geneSort_order = []
		geneSort_dict = {}
		for gene in gff_dict[chrom].keys():
			start = int(gff_dict[chrom][gene]["data"]["start"])
			same_start = sum([start == math.floor(i) for i in geneSort_order])
			if same_start > 0:
				start = start + same_start/100
			geneSort_order.append(start)
			geneSort_dict[str(start)] = gene

		geneSort_order.sort()
	
		for geneKey in geneSort_order:
			gene = geneSort_dict[str(geneKey)]
			gene = gff_dict[chrom][gene]
			printData(outfile, gene, False)

			mrna_ids = [i for i in gene.keys() if i != "data"]
			mrnaSort_order = []
			mrnaSort_dict = {}
			for mrna in mrna_ids:
				start = int(gene[mrna]["data"]["start"])
				same_start = sum([start == math.floor(i) for i in mrnaSort_order])
				if same_start > 0:
					start = start + same_start/100
				mrnaSort_order.append(start)
				mrnaSort_dict[str(start)] = mrna
		
			mrnaSort_order.sort()

			for mrnaKey in mrnaSort_order:
				mrna = mrnaSort_dict[str(mrnaKey)]
				mrna = gene[mrna]
				printData(outfile, mrna, True)

	outfile.close()


def tidyGFF(pre: str, gff: str, names: bool, out: str, splice: str, justify: int, no_chrom_id: bool, sort: bool, chrom_regex = False, contig_regex = False, source = None):

	# Guard against a path/directory being passed as the ID prefix (issue #10):
	# a '/' in `pre` lands inside every gene ID (e.g. "results/FILTER00000010"),
	# breaking downstream tools that use IDs as filenames/keys.
	if pre is None or "/" in str(pre):
		raise ValueError(f"TidyGFF prefix (pre) must be a non-path token, got: {pre!r}")

	if sort:
		sortGFF(loadGFF(gff), "preSort.tidyGFF.gff", chrom_regex, contig_regex)
	
	justify = justify - 1
	
	if not no_chrom_id:
		if sort:
			sorted = pd.read_csv("preSort.tidyGFF.gff", sep='\t', comment='#', header=None, names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
		else:
			sorted = pd.read_csv(gff, sep='\t', comment='#', header=None, names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
		# Get the total number of true chromosomes
		if chrom_regex:
			total_chroms = sorted.loc[sorted.seqid.str.contains(chrom_regex), "seqid"].unique()
		else:
			total_chroms = sorted.loc[sorted.seqid.str.contains(r'(^Chr)|(^chr)|(^LG)|(^Ch)|(^\d)'), "seqid"].unique()
		# Extract the numeric part of each chromosome name, skipping names with no
		# digits (ChrM/ChrC) or suffixes so int("") does not crash (issue #20.11).
		chrom_nums = [int(re.sub(r"\D", "", i)) for i in total_chroms if re.sub(r"\D", "", i)]
		total_chroms = max(chrom_nums) if chrom_nums else 1
		if total_chroms >= 10:
			chrom_just = 2
		else:
			chrom_just = 1

	if sort:
		read_file = open("preSort.tidyGFF.gff", 'r')
	else:
		read_file = open(gff, 'r')

	map_file = open(f"{os.path.basename(gff)}.map", 'w')
	out_file = open(out, 'w')

	chrom = None
	seq_count = 1
	transcript_count = 1
	feature_counts = {}   # per-transcript {feature_type: count}, reset at each mRNA
	seen_ids = set()      # global ID uniqueness guard (GFF3 requires unique feature IDs)

	def _register(full_id):
		"""Track emitted IDs and guarantee global uniqueness (issue #11).

		Child features are numbered per (transcript, type) so collisions cannot
		normally occur, but if one ever does (e.g. malformed input) append a
		numeric disambiguator instead of emitting a duplicate ID.
		"""
		if full_id not in seen_ids:
			seen_ids.add(full_id)
			return full_id
		n = 2
		while f"{full_id}_{n}" in seen_ids:
			n += 1
		dedup = f"{full_id}_{n}"
		print(f"WARNING: duplicate feature ID {full_id!r} -> {dedup!r}")
		seen_ids.add(dedup)
		return dedup

	if contig_regex:
		contig_regex = contig_regex.split(",")

	for line in read_file:
	
		line = line.strip()
		if not line or line.startswith("#"):
			continue

		line = line.split("\t")
		# 9-field normalisation, matching loadGFF (issue #20.8): rejoin a column-9
		# value that an embedded TAB split into extra fields; skip short lines.
		if len(line) < 9:
			continue
		if len(line) > 9:
			line = line[:8] + [" ".join(line[8:])]

		if not no_chrom_id:
			if chrom != line[0]:
				seq_count = 1
				chrom = line[0]
	
		if source:
			line[1] = source

		if names:
			line[8] = re.sub(r"Name=[a-zA-Z\d\.\-\_\%]*;{0,1}", "", line[8])

		if line[2] == "gene":
			if no_chrom_id:
				seq = f"{pre}{str(seq_count).zfill(justify)}0"
			else:
				new_chrom = upgradeChromosome(chrom, chrom_just, chrom_regex, contig_regex)
				seq = f"{pre}{new_chrom}{str(seq_count).zfill(justify)}0"
			seq_count += 1
			transcript_count = 1
			seq = _register(seq)
			line[8], old_id = replaceID(seq, line[8], "")
			map_file.write("\t".join(["gene", old_id, seq]) + "\n")
		elif line[2] == "mRNA":
			base = getID(line[8])
			line[8] = replaceParent(seq, line[8])
			transcript_seq = seq + f".{splice}{transcript_count}"
			transcript_count += 1
			feature_counts = {}   # restart child numbering for this transcript
			transcript_seq = _register(transcript_seq)
			line[8], old_id = replaceID(transcript_seq, line[8], "")
			map_file.write("\t".join(["mRNA", old_id, transcript_seq]) + "\n")
		else:
			line[8] = replaceParent(transcript_seq, line[8])
			# Number each child feature per (transcript, type) from the authoritative
			# type column, e.g. .cds1/.cds2/.exon1 (issue #11). The former getSuffix
			# path derived the suffix from the source ID and gave every CDS of a
			# multi-CDS mRNA the identical ".cds" (167k duplicate IDs in the wild),
			# and crashed with TypeError when the source ID lacked a cds/exon/utr marker.
			ftype = line[2].lower()
			feature_counts[ftype] = feature_counts.get(ftype, 0) + 1
			full_id = _register(f"{transcript_seq}.{ftype}{feature_counts[ftype]}")
			suffix = full_id[len(transcript_seq):]
			line[8], old_id = replaceID(transcript_seq, line[8], suffix)
			map_file.write("\t".join([line[2], old_id, full_id]) + "\n")

		out_file.write("\t".join(line) + "\n")

	out_file.close()	

if __name__ == "__main__":
	ap = argparse.ArgumentParser(
		prog='TidyGFF.py',
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=('Reformat and validate GFF files.')
	)
	ap.add_argument('pre', type=str, help="Prefix to be associated with gene IDs")
	ap.add_argument('gff', type=str, help="A GFF file")
	ap.add_argument('--out', type=str, default="tidy.gff", help="Output file name")
	ap.add_argument('--remove-names', action='store_true', help="Remove names from attributes")
	ap.add_argument('--splice-name', type=str, default="", help="Splice variant labels. For example, 'mRNA' will add [mRNA1, mRNA2, mRNA3 ...]")
	ap.add_argument('--justify', type=int, default=8, help="Number of digits in each gene ID")
	ap.add_argument('--no-chrom-id', action='store_true', help="If set, gene IDs will not be numbered based on chromosome")
	ap.add_argument('--sort', action='store_true', help="Perform sorting by chromosome and start coordinates")
	ap.add_argument('--chrom-regex', default=False, help=r"Provide regex for chromosome prefixes. Prefixes Chr, chr, LG, Ch, ^\d are automatically detected.")
	ap.add_argument('--contig-regex', default=False, help=r"Provide regex for contig/scaffold prefixes and the desired gene name prefix. e.g 'HiC_scaffold_(\d+$),Scaf' results in prefixScaf(\$d+)G000010")
	ap.add_argument('--source', type=str, default=None, help="Value for GFF column 2")
	args = ap.parse_args()

	pre = args.pre
	gff = args.gff
	names = args.remove_names
	out = args.out
	splice = args.splice_name
	justify = args.justify - 1
	no_chrom_id = args.no_chrom_id
	
	tidyGFF(args.pre, args.gff, args.remove_names, args.out, args.splice_name, args.justify , args.no_chrom_id, args.sort, args.chrom_regex, args.contig_regex, args.source)
	
	print(f"\nTidy GFF file written to {out}.")
