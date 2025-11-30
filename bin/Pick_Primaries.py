import pandas as pd
import argparse
import os
import re

parser = argparse.ArgumentParser(
			prog = 'Pick_Primaries',
			description = 'Pick out primary transcripts based on TransDecoder scores. Use ...transdecoder_dir/longest_orfs.cds.scores to pick primary transcripts.')

parser.add_argument("scores", type = str, help = "Transdecoder longestOrf.cds.scores file")
parser.add_argument("--gff", type = str, default = "", help = "Add transdecoder scores to the attributes of a gff3 file")
parser.add_argument("--cds", type = str, default = "", help = "Filter a fasta file of CDS sequences to provide primary transcripts")
args = parser.parse_args()

scores = pd.read_csv(args.scores, sep="\t", usecols=["#acc", "score_1"])
scores = scores[scores["#acc"].str.contains("\.p1$")].reset_index(drop=True)
scores.loc[:, "gene"] = scores.loc[:, "#acc"].str[:-6]
scores.loc[:, "transcript"] = scores.loc[:, "#acc"].str[:-3]

primary = scores[scores.groupby(['gene'])['score_1'].transform(max) == scores["score_1"]].drop_duplicates(["gene"])

primary["transcript"].to_csv("primary_transcripts.txt", sep = "\t", header=False, index=False)
print("Primary transcript ids written to: primary_transcripts.txt")

if args.cds != "":
	print("Filtering primary cds seqeunces ...")
	os.system(f"seqkit grep -n -r -f primary_transcripts.txt {args.cds} > primary_transcripts.fa")
	print("Primary transcript sequences written to: primary_transcripts.fa")

print("Updating GFF...")
if args.gff != "":
	out = open("primary_transcripts.gff3", "w")

	with open(args.gff, "r") as g:
		for line in g:
			line = line.strip()

			if line.startswith("#"):
				out.write(line + "\n")
				continue

			line = line.split("\t")

			if line[2] == "mRNA":
				tID = re.search("ID=([a-zA-Z\_\-\d\.]+);", line[8]).group(1)
				score = scores[scores["transcript"]==tID]["score_1"].values[0]
				line[8] = re.sub(";$", "", line[8]) + f";transdecoder_orf_score={score}"
				if tID in primary["transcript"].values:
					line[8] = line[8] + ";primary_transcript=True"

			out.write("\t".join(line) + "\n")
	
	out.close()
	print("Updated gff written to: primary_transcripts.gff3")
