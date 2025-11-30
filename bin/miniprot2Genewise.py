import re, sys, os, gzip
import pandas as pd

genome_file = sys.argv[1] # Genome file
protein_file = sys.argv[2] # Protein file
miniprot_file = sys.argv[3] # Miniprot file
output_dir = sys.argv[4] if len(sys.argv) > 4 else "results/GETA/homolog/geneRegion_genewise.tmp"  # Output directory

# Read fasta file into a dictionary (supports gzipped files)
def readFasta(path:str, sep=" " ,index=0) -> dict:
	seq = {}
	id = None
	lines = []

	# Handle gzipped files
	open_func = gzip.open if path.endswith('.gz') else open
	mode = 'rt' if path.endswith('.gz') else 'r'

	with open_func(path, mode) as file:	
		while True:
			line = file.readline()
			if not line:
				break

			line = line.strip()
			if line.startswith('>'):
				if lines:
					seq[id] = "".join(lines)
				id = line[1:]
				try:
					id = id.split(sep)[index].strip()
				except:
					id = id.split(" ")[0]
				lines = []
			else:
				lines.append(line)
	
	if id and lines:
			seq[id] = "".join(lines)

	return seq

# Read miniprot gff regions
miniprot_regions = {}
with open(miniprot_file, 'r') as file:
	for line in file:
		if line.startswith("#"):
			continue
		
		line = line.strip().split("\t")

		if line[2] == "mRNA":
			target = re.search("Target=(\S+)", line[8]).group(1)
			#if "|" in target:
			#	target = target.split("|")[1]
			miniprot_regions[target] = {
				"chr": line[0],
				"start": line[3],
				"end": line[4],
				"strand": line[6]
			}


# Read genome and protein_db into dictionary
genome = readFasta(genome_file)
proteins = readFasta(protein_file)

# Create directory if it doesn't exist
tmp_dir = output_dir
if not os.path.exists(tmp_dir):
	os.makedirs(tmp_dir)

# Process unique protein-gene mappings
fna = open(f"{tmp_dir}/group0.fna", 'w')
faa = open(f"{tmp_dir}/group0.faa", 'w')
for i, target in enumerate(miniprot_regions):
	seq_id = miniprot_regions[target]["chr"]
	start = miniprot_regions[target]["start"]
	end = miniprot_regions[target]["end"]
	strand = miniprot_regions[target]["strand"]
	if strand == "+":
		strand = "plus"
	else:
		strand = "minus"
	homolog_id = target

	if i % 100 == 0:
		fna.close()
		faa.close()
		command_group = i
		fna = open(f"{tmp_dir}/group{command_group}.fna", 'w')
		faa = open(f"{tmp_dir}/group{command_group}.faa", 'w')

	length = int(end) - int(start) + 1
	start = max(0, int(start) - 1)
	seq = genome[seq_id][start:start+length]

	fna.write(f">{seq_id}.{start}.{end}.{strand}\n{seq}\n")

	faa.write(f">{seq_id}.{start}.{end}.{strand}!{homolog_id}\n{proteins[homolog_id]}\n")

fna.close()
faa.close()
