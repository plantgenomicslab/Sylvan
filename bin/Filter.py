import os,re,sys,argparse
from pprint import pprint
import pandas as pd
import TidyGFF
from sklearn.ensemble import RandomForestClassifier
import numpy as np

'''
Usage: filter_test.py annotated.gff3 fpkm_cutoff coverage_cutoff augustus_cutoff helixer_cutoff repeat_cutoff blast_pident blast_qcovs rex_pident blast_qcovs

Script filters gene models with poor evidence from a
gff3 file and creates a filtered gff3. It relies on six previously 
generated files:
	- FILTER/BLASTP.OUT.TMP
	- FILTER/BLASTP.OUT.Rex
	- FILTER/cov.bed
	- FILTER/pfam.out
	- FILTER/RSEM.genes.results
	- FILTER/augustus_coverage.bed
	- FILTER/helixer_coverage.bed
	- FILTER/repeat_coverage.bed
'''

def read_rsem_file(file_path):
	return pd.read_csv(file_path, delimiter="\t")

def read_blast_file(file_path):
	return pd.read_csv(file_path, header=None, sep="\t", usecols=[0,2,12], names=["transcript_id","pident", "qcovs"])

def read_cov_file(file_path):
	return pd.read_csv(file_path, header=None, sep="\t", usecols=[0,6])

def read_pfam_file(file_path):
	return pd.read_csv(file_path, header=None, delim_whitespace=True, comment='#', usecols=[0])

def read_gff_file(file_path):
	return pd.read_csv(file_path, header=None, sep="\t")

def read_abInitio_cov(file_path):
	return pd.read_csv(file_path, header=None, sep="\t", usecols=[4,8])

def filter_genes(tpm_cutoff, cov_cutoff, blast_pident, blast_qcovs, rex_pident, rex_qcovs, output_dir="FILTER"):
	# Collect evidence and perform cutoffs
	## RSEM results
	rsem_file = os.path.join(output_dir, "rsem_outdir", "RSEM.isoforms.results")
	rsem_data = read_rsem_file(rsem_file)
	filter = rsem_data.loc[:, ("transcript_id", "TPM")]
	data = rsem_data.loc[:, ("transcript_id", "TPM")]
	filter.loc[:,"TPM"] = rsem_data["TPM"] > float(tpm_cutoff)

	## BlastP results
	blast_file = os.path.join(output_dir, "BLASTP.OUT.TMP")
	blast_data = read_blast_file(blast_file)
	blast_data.drop_duplicates(subset=["transcript_id"], keep='first', inplace=True, ignore_index=True)
	blast_data["blast_pident"] = blast_data["pident"]
	blast_data["blast_qcovs"] = blast_data["qcovs"]
	data = data.merge(blast_data.loc[:, ("transcript_id", "blast_pident", "blast_qcovs")], on="transcript_id", how="left")
	blast_data["blast_pident"] = blast_data["pident"] > float(blast_pident)*100
	blast_data["blast_qcovs"] = blast_data["qcovs"] > float(blast_qcovs)*100 
	filter = filter.merge(blast_data.loc[:, ("transcript_id", "blast_pident", "blast_qcovs")], on="transcript_id", how="left")
	filter.fillna(False, inplace=True)
	data.fillna(0, inplace=True)

	## Coverage results
	cov_file = os.path.join(output_dir, "cov.bed")
	cov_data = read_cov_file(cov_file)
	cov_data.columns = ["transcript_id", "COVERAGE"]
	data = data.merge(cov_data, on="transcript_id", how="left")
	
	cov_data["COVERAGE"] = cov_data["COVERAGE"] > float(cov_cutoff)
	filter = filter.merge(cov_data, on="transcript_id", how="left")
	
	## Pfam results
	pfam_file = os.path.join(output_dir, "pfam.out")
	pfam_data = read_pfam_file(pfam_file)
	pfam_data = pfam_data[0].unique()
	pfam_data = pd.DataFrame(pfam_data, columns=["transcript_id"])
	pfam_data["PFAM"] = True
	data = data.merge(pfam_data, on="transcript_id", how="left")
	filter = filter.merge(pfam_data, on="transcript_id", how="left")
	filter.drop_duplicates(subset=["transcript_id"], keep='first', inplace=True, ignore_index=True)
	filter.fillna(False, inplace=True)
	data.fillna(False, inplace=True)

	## Ab Initio coverage results
	augustus_cov = read_abInitio_cov(os.path.join(output_dir, "augustus_coverage.bed"))
	helixer_cov = read_abInitio_cov(os.path.join(output_dir, "helixer_coverage.bed"))
	augustus_cov.columns = ["transcript_id", "AUGUSTUS"]
	helixer_cov.columns = ["transcript_id", "HELIXER"]
	
	data = data.merge(augustus_cov, on="transcript_id", how="left")
	data = data.merge(helixer_cov, on="transcript_id", how="left")
	
	augustus_cov["AUGUSTUS"] = augustus_cov["AUGUSTUS"] > float(augustus_cutoff)
	helixer_cov["HELIXER"] = helixer_cov["HELIXER"] > float(helixer_cutoff)
	filter = filter.merge(augustus_cov, on="transcript_id", how="left")
	filter = filter.merge(helixer_cov, on="transcript_id", how="left")
	filter.fillna(False, inplace=True)
	data.fillna(0, inplace=True)

	## Repeat coverage results
	repeat_cov = read_abInitio_cov(os.path.join(output_dir, "repeat_coverage.bed"))
	repeat_cov.columns = ["transcript_id", "REPEAT"]
	data = data.merge(repeat_cov, on="transcript_id", how="left")
	
	repeat_cov["REPEAT"] = repeat_cov["REPEAT"] < float(repeat_cutoff)
	filter = filter.merge(repeat_cov, on="transcript_id", how="left")
	filter.fillna(False, inplace=True)
	data.fillna(0, inplace=True)

	## Single Exons
	singleExons = TidyGFF.singleExonGenes(TidyGFF.loadGFF(sys.argv[1]))
	singleExons["singleExon"] = True
	data = data.merge(singleExons, on="transcript_id", how="left")
	filter = filter.merge(singleExons, on="transcript_id", how="left")
	filter.fillna(False, inplace=True)
	data.fillna(False, inplace=True)

	## RexDB Blastp
	rexdb = os.path.join(output_dir, "BLASTP.OUT.Rex")
	rexdb = read_blast_file(rexdb)
	rexdb.drop_duplicates(subset=["transcript_id"], keep='first', inplace=True, ignore_index=True)
	rexdb["rex_pident"] = rexdb["pident"]
	rexdb["rex_qcovs"] = rexdb["qcovs"]
	data = data.merge(rexdb.loc[:, ("transcript_id", "rex_pident", "rex_qcovs")], on="transcript_id", how="left")
	rexdb["rex_pident"] = rexdb["pident"] > float(rex_pident)*100
	rexdb["rex_qcovs"] = rexdb["qcovs"] > float(rex_qcovs)*100
	filter = filter.merge(rexdb.loc[:, ("transcript_id", "rex_pident", "rex_qcovs")], on="transcript_id", how="left")
	filter.fillna(False, inplace=True)
	data.fillna(0, inplace=True)
	
	## LncRNA prediction
	lncrna = pd.read_csv(os.path.join(output_dir, "lncrna_predict.csv"), comment='#', usecols=[0,25], names=["transcript_id","LncRNA_predict"])
	lncrna["transcript_id"] = lncrna["transcript_id"].str.split(" ").str.get(0)
	lncrna["LncRNA_predict"] = lncrna["LncRNA_predict"] == "lncrna"
	data = data.merge(lncrna, on="transcript_id", how="left")
	filter = filter.merge(lncrna, on="transcript_id", how="left")
	filter.fillna(False, inplace=True)
	data.fillna(False, inplace=True)

	# Identify strong keep hits while maintaining as much feature diversity as possible
	data["keep"] = pd.Series(((filter["COVERAGE"] == True) & ((filter["blast_pident"] == True) & (filter["blast_qcovs"] == True) & (filter["PFAM"] == True) | ((filter["AUGUSTUS"] == True) & (filter["HELIXER"] == True)))) |
			  ((filter["singleExon"] == False) & (filter["LncRNA_predict"] == False) & (filter["COVERAGE"] == True) & (filter["TPM"] == True)) |
			  ((filter["singleExon"] == False) & (filter["PFAM"] == True) & (filter["blast_pident"] == True) & (filter["blast_qcovs"] == True)) |
			  ((filter["AUGUSTUS"] == True) & (filter["HELIXER"] == True) & ((filter["rex_pident"] == False) & (filter["rex_qcovs"] == False)) & (filter["REPEAT"] == False) & (filter["LncRNA_predict"] == False) & (filter["singleExon"] == False)))

	# Identify strong discard hits while maintaining as much feature diversity as possible
	data["discard"] = pd.Series((filter["COVERAGE"] == False)  & (filter["AUGUSTUS"] == False) & (filter["HELIXER"] == False) & (filter["blast_pident"] == False) & (filter["blast_qcovs"] == False) & (filter["PFAM"] == False) |
		((filter["singleExon"] == True) & ((filter["blast_pident"] == False) & (filter["blast_qcovs"] == False)) & (filter["PFAM"] == False) & (filter["COVERAGE"] == False)) |
		(((filter["rex_pident"] == True) & (filter["rex_qcovs"] == True)) & (filter["REPEAT"] == True)) & ((filter["blast_pident"] == False) & (filter["blast_qcovs"] == False)) |
		(filter["LncRNA_predict"] == True) & (filter["singleExon"] == True) & (filter["PFAM"] == "False") & (filter["AUGUSTUS"] == False) & (filter["HELIXER"] == False))
	
	data["label"] = data.apply(lambda x: "Discard" if x["discard"] else ("Keep" if x["keep"] else "None"), axis=1)
	
	return data.drop("discard", axis=1).drop("keep", axis=1)

def filter_gff(gff_data, keep):
	gff_data["transcript_id"] = gff_data[8].apply(getGeneId)
	names = keep['New_ID']
	names = pd.concat([names, names.str.split(".").str.get(0)])
	to_keep = gff_data["transcript_id"].isin(names)
	gff_keep = gff_data[to_keep]
	gff_discard = gff_data[~to_keep]
	return gff_keep, gff_discard

def getGeneId(attr):
	m = re.search("ID=(FILTER\d+)(\.t[0-9]+){0,1}(\.[a-zA-Z\d\_\.]+){0,1};", attr)
	match = m.group(1)
	try:
		return(match + m.group(2))
	except:
		return(m.group(1))

def getParent(attr):
	m = re.search("Parent=[a-zA-Z\d\._-]*;", attr)
	if m:
		m = m.group(0).replace("Parent=", "").replace(";", "")
	else:
		m = None
	return(m)

def checkKeep(id, names):
	val = sum(names.str.match(id)) >= 1
	return(val)

def formatKeepIDs(id):
	base = re.search("^(FILTER\d+)\.*", id).group(1)
	return(base)

def getGene(data, starts, name_vec):
	ns = data.loc[(data[2]=='mRNA') & (data[8].str.startswith(starts))]
	ns = ns.loc[ns['transcript_id'].isin(name_vec)]
	ns = ns.iloc[:,8].apply(getParent)
	return(ns)

def semiSupRandomForest(data, predictors, busco_table, num_trees, seed=None, recycle_prob=0.95, maxiter=5):
	# data        : initialized data with labeled and unlabled cases
	# predictors  : Number of predictors to consider when splitting each node
	# busco_table : Busco table from gff being filtered
	# num_trees   : number of trees in the forest
	# seed        : random seed for reproducibility
	# recycle_prob: The prediction probability required to label an un-labeled gene model

	busco = pd.read_csv(busco_table, sep="\t", comment="#", usecols=(0,1,2), names=["busco_id","status","transcript_id"])
	busco = busco.loc[busco["status"] != "Missing"].reset_index(drop=True)
	
	train = data.loc[data["label"] != "None"].set_index("transcript_id")
	test = data.loc[data["label"] == "None"].set_index("transcript_id")

	#for i in range(0, 10):
	x_train = train.loc[:, train.columns != "label"]
	y_train = train.loc[:, train.columns == "label"]
	x_test = test.loc[:, test.columns != "label"]
	
	stop = False
	iter = 0
	process = {
		"kept":[],
		"discarded":[],
		"kept_buscos":[],
		"discarded_buscos":[],
		"OOB":[]
	}
	while stop == False:
		if iter >= maxiter:
			print(f"Training stopped at a maximum of {maxiter} iterations.")
			break

		discarded = y_train[y_train["label"] == "Discard"].index
		kept = y_train[y_train["label"] == "Keep"].index
		
		process["discarded"].append(len(discarded))
		process["kept"].append(len(kept))
		
		kept_buscos = []
		dis_buscos = []
		for i in range(len(busco)): # TODO: also report total kept buscos instead of only the buscos labled as 'keep'
			if busco.loc[i, "transcript_id"] in kept and busco.loc[i, "busco_id"] not in kept_buscos:
				kept_buscos.append(busco.loc[i, "busco_id"])
			elif busco.loc[i, "transcript_id"] in discarded and busco.loc[i, "busco_id"] not in dis_buscos:
				dis_buscos.append(busco.loc[i, "busco_id"])
		
		process["discarded_buscos"].append(len(dis_buscos))
		process["kept_buscos"].append(len(kept_buscos))

		model = RandomForestClassifier(n_estimators = num_trees, 
				 max_features = predictors, 
				 random_state = seed,
				 oob_score=True)
		model.fit(x_train, y_train.to_numpy().reshape(-1,))
		process["OOB"].append(1-model.oob_score_)
		
		yhat = pd.DataFrame(model.predict(x_test), columns=["prediction"])
		yprob = pd.DataFrame(model.predict_proba(x_test), columns=model.classes_)
		yhat = pd.concat([yhat, yprob], axis=1)
		yhat["pred_prob"] = yhat.apply(lambda x: x[x["prediction"]], axis=1)
		yhat["transcript_id"] = x_test.index
		
		recycle = yhat.loc[yhat["pred_prob"] > recycle_prob, "transcript_id"]
		if len(recycle) > 0:
			recycle = pd.DataFrame(recycle).merge(data, on="transcript_id", how="left")
			recycle["label"] = yhat["prediction"]
			recycle = recycle.set_index("transcript_id")
		
			x_train = pd.concat([x_train, recycle.loc[:, recycle.columns != "label"]], axis=0)
			y_train = pd.concat([y_train, recycle.loc[:, recycle.columns == "label"]], axis = 0)
			x_test  = x_test.loc[[x_test.index[i] not in recycle.index for i in range(len(x_test))]]
			iter += 1
		else:
			stop = True
			break
	
	return pd.DataFrame(kept, columns=["transcript_id"]), process


if __name__ == "__main__":
	# Parse arguments
	parser = argparse.ArgumentParser(
			prog = 'Filter.py',
			description = 'Filter gene models using a semi-supervised random forest.')

	parser.add_argument("gff", type = str, help = "GFF file to be filtered")
	parser.add_argument("busco", type = str, help = "BUSCO output for monitoring filter")
	parser.add_argument("--tpm", type = float, default = 3.0, help = "TPM threshold for initial filtration (Default: 3)")
	parser.add_argument("--rcov", type = float, default = 0.5, help = "RNAseq coverage threshold for initial filtration (Default: 0.5)")
	parser.add_argument("--acov", type = float, default = 0.8, help = "Augustus coverage threshold for initial filtration (Default: 0.8)")
	parser.add_argument("--hcov", type = float, default = 0.8, help = "Helixer coverage threshold for initial filtration (Default: 0.8)")
	parser.add_argument("--repcov", type = float, default = 0.5, help = "Repeat coverage threshold for initial filtration (Default: 0.5)")
	parser.add_argument("--blast-pident", type = float, default = 0.6, help = "BLASTp percent identity threshold for initial filtration (Default: 0.6)")
	parser.add_argument("--blast-qcovs", type = float, default = 0.6, help = "BLASTp query coverage threshold for initial filtration (Default: 0.6)")
	parser.add_argument("--rex-pident", type = float, default = 0.6, help = "REXdb percent identity threshold for initial filtration (Default: 0.6)")
	parser.add_argument("--rex-qcovs", type = float, default = 0.6, help = "REXdb query coverage threshold for initial filtration (Default: 0.6)")
	parser.add_argument("--predictors", type = int, default = 6, help = "Number of predictors used for tree splitting during random forest tree building (Default: 6)")
	parser.add_argument("--recycle", type = float, default = 0.95, help = "Predicted accuracy required for an observation to be added to the model for successive iterations of the random forest (Default: 0.95)")
	parser.add_argument("--trees", type = int, default = 100, help = "Number of trees in the random forest (Default: 100)")
	parser.add_argument("--max-iter", type = int, default = 5, help = "Maximum number of random forest re-training iterations (Default: 5)")
	parser.add_argument("--seed", type = int, default = 123, help = "Random seed for reproducibility (Default: 123)")
	parser.add_argument("--chromRegex", type = str, default = "", help = "Regular expression to match chromosome prefixes (Default: '')")
	parser.add_argument("--output", type = str, default = "filter.gff3", help = "Output file path (Default: filter.gff3)")
	parser.add_argument("--output-dir", type = str, default = "FILTER", help = "Output directory for intermediate files (Default: FILTER)")
	args = parser.parse_args()

	gff = args.gff
	busco = args.busco
	tpm_cutoff = args.tpm
	cov_cutoff = args.rcov
	augustus_cutoff = args.acov
	helixer_cutoff = args.hcov
	repeat_cutoff = args.repcov
	blast_pident = args.blast_pident
	blast_qcovs = args.blast_qcovs
	rex_pident = args.rex_pident
	rex_qcovs = args.rex_qcovs

	output_dir = args.output_dir
	output_file = args.output

	# Apply initial round of labeling
	data = filter_genes(tpm_cutoff, cov_cutoff, blast_pident, blast_qcovs, rex_pident, rex_qcovs, output_dir)
	data.to_csv(f"{output_dir}/data.tsv", sep = "\t", index=False)
	# Apply semi-supervised learning to classify gene models
	keep, report = semiSupRandomForest(data, args.predictors, busco, args.trees, seed=args.seed, recycle_prob=args.recycle, maxiter=args.max_iter)
	pprint(report)

	# Reformat GFF with easier IDs
	TidyGFF.tidyGFF(output_dir, gff, True, "tidy.gff", "t", 8 , True, True, chrom_regex=args.chromRegex)
	map_file = pd.read_csv(f"{os.path.basename(gff)}.map", sep = "\t", names = ["transcript_id", "New_ID"], usecols = [1,2])
	gff_data = read_gff_file("tidy.gff")


	keep = keep.merge(map_file, on="transcript_id", how="left")
	keep["basenames"] = keep['New_ID'].apply(formatKeepIDs)

	# Filter gff
	gff_keep, gff_discard = filter_gff(gff_data, keep)

	# Write output
	gff_keep.iloc[:, 0:9].to_csv(output_file, sep = "\t", index=False, header=False)
	gff_discard.iloc[:, 0:9].to_csv(f"{output_dir}/discard.gff3", sep = "\t")

	# Display stats
	print(f"\nTotal genes	: {sum(gff_data[2] == 'gene')}")
	print(f"Total mRNA	: {sum(gff_data[2] == 'mRNA')}\n")
	print(f"\nGenes kept	 : {sum(gff_keep[2] == 'gene')}")
	print(f"mRNA kept	 : {sum(gff_keep[2] == 'mRNA')}\n")
	print(f"\nGenes discarded: {sum(gff_discard[2] == 'gene')}")
	print(f"mRNA discarded: {sum(gff_discard[2] == 'mRNA')}\n")

	print(f"Filtered gene models written to: {output_file}")
	print(f"Discarded gene models written to: {output_dir}/discard.gff3\n")

	print(f"Filter data written to: {output_dir}/discard_data.tsv   {output_dir}/keep_data.tsv")
	print(f"See {gff}.map for associations with original IDs.")
