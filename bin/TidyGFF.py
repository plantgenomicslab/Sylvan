import argparse, re, math, os, subprocess
import pandas as pd

def getParent(attr: str) -> str:
	return re.search("Parent=([a-zA-Z\d\\.\-_:]*);{0,1}", attr).group(1)

def getID(attr:str) -> str:
	return re.search("ID=([a-zA-Z\d\\.\-_:]*);{0,1}", attr).group(1)

def findChroms(chrom: str, chrom_regex = False) -> str:
	if chrom_regex:
		if re.search(chrom_regex, chrom):
			pre = re.search(chrom_regex, chrom).group(0)
	elif re.search("(^Chr)|(^chr)|(^LG)|(^Ch)|(^\d)", chrom):
		pre = re.search("(^Chr)|(^chr)|(^LG)|(^Ch)|(^\d)", chrom).group(0)
	return re.sub(pre, "", chrom)

def replaceParent(id: str, attr: str) -> str:
	replace = f"Parent={id};"
	return re.sub("Parent=[a-zA-Z\d\.\-\_]*;{0,1}", replace, attr)

def replaceID(id: str, attr: str, suf: str) -> str:
	replace = f"ID={id}{suf};"
	replace = re.sub("ID=[a-zA-Z\d\.\-\_]*;{0,1}", replace, attr)
	old = re.search("ID=([a-zA-Z\d\.\-\_]*);{0,1}", attr).group(1)
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
		chrom_num = int(chrom.replace(m[0], ""))
		new_name = str(chrom_num).zfill(just) + "G"
		return(new_name)
	elif re.search("^\d", chrom):
		new_name = chrom.zfill(just) + "G"
		return(new_name)
	
	if contig_regex:
		if re.search(contig_regex[0], chrom):
			contig_num = re.search("[\d\_\-]+", chrom)[0]
			new_name = contig_regex[1] + str(contig_num) + "G"
			return(new_name)
	
	return(chrom)

def getSuffix(id: str) -> str:
	id = id.lower()
	suf = ["cds", "exon", "intron", "utr5p", "utr3p", "five_prime_utr", "three_prime_utr"]
	for suffix in suf:
		if re.search(suffix, id):
			return(re.search(f"{suffix}\.{{0,1}}\d*", id)[0].lower())

def loadGFF(gff: str) -> dict:
	with open(gff, 'r') as infile:
		gff_dict = {}
		missing_parent = False

		for line in infile:
			if line.startswith("#") or not line.strip():
				continue

			line = line.strip().split("\t")

			#### For loading mikado output ###
			if (line[2] == "superlocus"):
				continue
			
			if (line[2] == "ncRNA_gene"):
				protein_coding = False
			##################################
			
			chrom = line[0]
			id = getID(line[8])
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
				cur_gene = id
				gff_dict[chrom][id] = {"data":feature_data}
				protein_coding = True # Support for mikado
			elif line[2] == "mRNA":
				cur_mrna = id
				gene = getParent(line[8])
				missing_parent = False
		
				if gene not in gff_dict[chrom].keys():
					print(f"WARNING: mRNA with id {id} missing parent {gene}")
					gff_dict[chrom][id] = {"data":feature_data}
					missing_parent = True
				else:
					if gene != cur_gene:
						print(f"WARNING: Could not parse parents for {line[2]} type feature with id {id}")
				
					gff_dict[chrom][gene][id] = {"data":feature_data}
			else:
				if not protein_coding: #Support for mikado
					continue
				mrna = getParent(line[8])
				if missing_parent:
					gff_dict[chrom][mrna][id] = {"data":feature_data}
				else:
					gene = getParent(gff_dict[chrom][gene][mrna]["data"]["attributes"])
					if (mrna != cur_mrna) | (gene != cur_gene):
						print(f"WARNING: Could not parse parents for {line[2]} type feature with id {id}")

					gff_dict[chrom][gene][mrna][id] = {"data":feature_data}
	return(gff_dict)

def singleExonGenes(gff_dict):
	singleExons = []
	for chrom in gff_dict:
		for gene in gff_dict[chrom]:
			mrna_ids = [i for i in gff_dict[chrom][gene].keys() if i != "data"]
			for mrna in mrna_ids:
				cds = 0
				feature_ids = [i for i in gff_dict[chrom][gene][mrna].keys() if i != "data"]
				for feature in feature_ids:
					if gff_dict[chrom][gene][mrna][feature]["data"]["type"] == "CDS":
						cds +=1
				if cds <= 1:
					singleExons.append(mrna)
	return(pd.DataFrame(singleExons, columns=["transcript_id"]))
			
def sortGFF(gff_dict: dict, output: str, chrom_regex = False, contig_regex = False) -> None:
	outfile = open(output, 'w')

	if not chrom_regex:
		chrom_regex = "(^Chr)|(^chr)|(^LG)|(^Ch)|(^\d)"
	if not contig_regex:
		contig_regex = "(?! (^Chr)|(^chr)|(^LG)|(^Ch)|(^\d))"
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
			total_chroms = sorted.loc[sorted.seqid.str.contains('(^Chr)|(^chr)|(^LG)|(^Ch)|(^\d)'), "seqid"].unique()
		total_chroms = max([int(re.sub("[A-Za-z]+", "", i)) for i in total_chroms])
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

	if contig_regex:
		contig_regex = contig_regex.split(",")

	for line in read_file:
	
		line = line.strip()
		if not line or line.startswith("#"):
			continue
	
		line = line.split("\t")
		
		if not no_chrom_id:
			if chrom != line[0]:
				seq_count = 1
				chrom = line[0]
	
		if source:
			line[1] = source

		if names:
			line[8] = re.sub("Name=[a-zA-Z\d\.\-\_\%]*;{0,1}", "", line[8])

		if line[2] == "gene":
			if no_chrom_id:
				seq = f"{pre}{str(seq_count).zfill(justify)}0"
			else:
				new_chrom = upgradeChromosome(chrom, chrom_just, chrom_regex, contig_regex)
				seq = f"{pre}{new_chrom}{str(seq_count).zfill(justify)}0"
			seq_count += 1
			transcript_count = 1
			line[8], old_id = replaceID(seq, line[8], "")
			map_file.write("\t".join(["gene", old_id, seq]) + "\n")
		elif line[2] == "mRNA":
			base = getID(line[8])
			line[8] = replaceParent(seq, line[8])
			transcript_seq = seq + f".{splice}{transcript_count}"
			transcript_count += 1
			line[8], old_id = replaceID(transcript_seq, line[8], "")
			map_file.write("\t".join(["mRNA", old_id, transcript_seq]) + "\n")
		else:
			line[8] = replaceParent(transcript_seq, line[8])
			id = getID(line[8])
			suffix = "." + getSuffix(id)
			line[8], old_id= replaceID(transcript_seq, line[8], suffix)
			map_file.write("\t".join([line[2], old_id, transcript_seq + suffix]) + "\n")

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
	ap.add_argument('--chrom-regex', default=False, help="Provide regex for chromosome prefixes. Prefixes Chr, chr, LG, Ch, ^\d are automatically detected.")
	ap.add_argument('--contig-regex', default=False, help="Provide regex for contig/scaffold prefixes and the desired gene name prefix. e.g 'HiC_scaffold_(\d+$),Scaf' results in prefixScaf(\$d+)G000010")
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
