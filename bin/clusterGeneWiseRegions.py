import duckdb
import sys
import os
import pandas as pd
from fasta_utils import readFasta

genome_file = sys.argv[1] # Genome file
protein_file = sys.argv[2] # Protein file
gff_file = sys.argv[3] # Helixer gff
region_file = sys.argv[4] # GETA/homolog/homolog_gene_region.tab
segment_threshold = sys.argv[5] # Segment threshold

# Parse gene:transcript relationships
def parse_gff(gff:str) -> dict:
	genes = {}
	with open(gff, 'r') as file:
		for line in file:
			if line.startswith("#") or line == "\n":
				continue
			line = line.strip().split("\t")
			seq_id = line[8].split(";")[0].split("=")[1]
			if line[2] == "gene":
				gene_id  = seq_id
			if line[2] == "mRNA":
				genes[seq_id] = gene_id
	return genes

# Read in gene regions
regions = duckdb.read_csv(region_file)

# Parse gene:transcript relationships for neighbor species proteins
gene_map = parse_gff(gff_file)
gene_map = pd.DataFrame.from_dict([gene_map]).T.reset_index()
gene_map.columns = ["transcript_id", "gene_id"]

duckdb.sql("CREATE SEQUENCE clustered_regions_seq START 1")

# Collapse protein mappings from the same genes that are within the segment threshold
clustered_regions = duckdb.sql(
"WITH regions_cte AS (" +
	"SELECT column0 as chromosome, column1 as start, column2 as fin, column3 as prot, column4 as strand " +
	"FROM regions), " +
"added_gene_id AS (" +
	"SELECT " +
		"r.chromosome, " +
		"r.prot, " +
		"r.start, " +
		"r.fin, " +
		"r.strand, " +
		"CASE WHEN gm.gene_id IS NULL THEN r.prot ELSE gm.gene_id END AS gene " +
	"FROM regions_cte AS r " +
	"LEFT JOIN gene_map AS gm " +
	"ON r.prot = gm.transcript_id), "
"GroupedProteins AS ( " +
	"SELECT " + 
		"r1.chromosome, " +
		"r1.prot AS prot, " +
		"r1.gene AS gene, " +
		"r1.strand AS strand, " +
		"r1.start AS start1, " + 
		"r1.fin AS fin1, " +  
		"r2.start AS start2, " + 
		"r2.fin AS fin2, " + 
	"FROM added_gene_id r1 " + 
	"LEFT JOIN added_gene_id r2 " + 
	"ON r1.chromosome = r2.chromosome " + 
		"AND r1.gene = r2.gene " + 
		"AND ((r1.start != r2.start) OR (r1.fin != r2.fin)) " +
		"AND (" + 
			f"(ABS(r1.start - r2.start) <= {segment_threshold}) " + 
			f"OR (ABS(r1.start - r2.fin) <= {segment_threshold}) " + 
			f"OR (ABS(r1.fin - r2.start) <= {segment_threshold}) " + 
			f"OR (ABS(r1.fin - r2.fin) <= {segment_threshold})) " + 
	"ORDER BY r1.chromosome, r1.gene, r1.start, r2.start) " 
"SELECT " +
	"chromosome, " +
	"gene, " +
	"prot, " +
	"strand, " +
	"start1, " +
	"fin1, " +
	"start2, " +
	"fin2, " +
	"LAG(fin1) OVER () AS prev_end_value " +
"FROM GroupedProteins ").to_df()

# Set cluster membership within genes based on segmentation threshold
clustered_regions["cluster"] = None
cluster = 0
prev_gene = None
for i, row in clustered_regions.iterrows():
	gene = row["gene"]
	
	if gene != prev_gene:
		new_cluster = True
	elif abs(row["start1"] - row["prev_end_value"]) > int(segment_threshold):
		new_cluster = True
	else:
		new_cluster = False
	
	if new_cluster:
		cluster += 1

	clustered_regions.loc[i, "cluster"] = cluster
	prev_gene = gene

# Aggregate the cluster members
aggregate_clusters = duckdb.sql(
"WITH AggregatedGroups AS ( " + 
	"SELECT " + 
		"chromosome, " + 
		"gene, " +
		"STRING_AGG(DISTINCT prot, ', ') AS grouped_proteins, " + 
		"MIN(LEAST(start1, start2)) AS group_start, " + 
		"MAX(GREATEST(fin1, fin2)) AS group_end, " +
		"strand " +
	"FROM clustered_regions " + 
	"GROUP BY chromosome, gene, strand, cluster) " + 
"SELECT " + 
    "chromosome, " + 
	"gene, " +
    "grouped_proteins, " + 
    "group_start, " + 
    "group_end, " +
	"strand " +
"FROM AggregatedGroups " + 
"ORDER BY chromosome, group_start; "
).to_df()

# Set cluster membership for overlapping genes
aggregate_clusters["cluster"] = None
cluster = 0
prev_group_end = 0
for i, row in aggregate_clusters.iterrows():
	if i > 0:
		prev_group_end = aggregate_clusters.loc[i-1, "group_end"]
	else:
		prev_group_end = 0

	if (row["group_start"] - prev_group_end) > 0:
		cluster += 1

	aggregate_clusters.loc[i, "cluster"] = cluster

# Collapse overlapping mappings from different genes
overlapped_regions = duckdb.sql(
"SELECT " +
	"chromosome, " +
	"cluster, " +
	"STRING_AGG(DISTINCT gene, ', ') AS grouped_genes, " +
	"STRING_AGG(DISTINCT grouped_proteins, ', ') AS grouped_proteins, " +
	"MIN(group_start) AS group_start, " +
	"MAX(group_end) AS group_end, " +
	"strand " +
"FROM aggregate_clusters " +
"GROUP BY chromosome, strand, cluster " +
"ORDER BY chromosome, group_start; "
).to_df()

# Read genome and protein_db into dictionary
genome = readFasta(genome_file)
proteins = readFasta(protein_file, sep="|", index=1)

# Create directory if it doesn't exist
tmp_dir = "GETA/homolog/geneRegion_genewise.tmp"
if not os.path.exists(tmp_dir):
	os.makedirs(tmp_dir)

# Process unique protein-gene mappings
fna = open(f"{tmp_dir}/group0.fna", 'w')
faa = open(f"{tmp_dir}/group0.faa", 'w')
for i, row in overlapped_regions.iterrows():
	seq_id = row["chromosome"]
	start = row["group_start"]
	end = row["group_end"]
	strand = row["strand"]
	homolog_ids = row["grouped_proteins"].split(", ")

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

	for homolog_id in homolog_ids:
		homolog_seq = proteins[homolog_id]
		faa.write(f">{seq_id}.{start}.{end}.{strand}!{homolog_id}\n{homolog_seq}\n")

fna.close()
faa.close()