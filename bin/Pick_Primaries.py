import pandas as pd
import argparse
import os
import re
import subprocess

parser = argparse.ArgumentParser(
			prog = 'Pick_Primaries',
			description = 'Pick out primary transcripts based on TransDecoder scores. Use ...transdecoder_dir/longest_orfs.cds.scores to pick primary transcripts.')

parser.add_argument("scores", type = str, help = "Transdecoder longestOrf.cds.scores file")
parser.add_argument("--gff", type = str, default = "", help = "Add transdecoder scores to the attributes of a gff3 file")
parser.add_argument("--cds", type = str, default = "", help = "Filter a fasta file of CDS sequences to provide primary transcripts")
args = parser.parse_args()

scores = pd.read_csv(args.scores, sep="\t", usecols=["#acc", "score_1"])
scores = scores[scores["#acc"].str.contains(r"\.p1$")].reset_index(drop=True)
# acc looks like "<gene>.t<N>.p1". The primary ORF suffix ".p1" is always 3
# chars, so transcript = acc[:-3]. The gene is transcript with its ".t<N>"
# suffix stripped -- do it with a regex, NOT a fixed str[:-6], which only
# removes ".t1.p1" and mis-splits multi-digit transcripts (".t10.p1" left a
# stray ".t1", splitting one gene into two "primaries").
scores.loc[:, "transcript"] = scores.loc[:, "#acc"].str[:-3]
scores.loc[:, "gene"] = scores.loc[:, "transcript"].str.replace(r"\.t\d+$", "", regex=True)

primary = scores[scores.groupby(['gene'])['score_1'].transform(max) == scores["score_1"]].drop_duplicates(["gene"])

primary["transcript"].to_csv("primary_transcripts.txt", sep = "\t", header=False, index=False)
print("Primary transcript ids written to: primary_transcripts.txt")

if args.cds != "":
	print("Filtering primary cds seqeunces ...")
	subprocess.run(f"seqkit grep -n -r -f primary_transcripts.txt {args.cds} > primary_transcripts.fa", shell=True, check=True)
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
				m = re.search(r"ID=([a-zA-Z\_\-\d\.]+);", line[8])
				if not m:
					out.write("\t".join(line) + "\n")
					continue
				tID = m.group(1)
				matched = scores.loc[scores["transcript"] == tID, "score_1"].values
				if len(matched) == 0:
					# mRNA with no TransDecoder score (e.g. no primary ORF):
					# leave it untouched instead of crashing on values[0].
					out.write("\t".join(line) + "\n")
					continue
				score = matched[0]
				line[8] = re.sub(";$", "", line[8]) + f";transdecoder_orf_score={score}"
				if tID in primary["transcript"].values:
					line[8] = line[8] + ";primary_transcript=True"

			out.write("\t".join(line) + "\n")
	
	out.close()
	print("Updated gff written to: primary_transcripts.gff3")
