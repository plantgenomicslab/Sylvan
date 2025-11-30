import pickle, sys, os, re, argparse, subprocess
import duckdb
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

class genome:
	def __init__(self, db="CombineGFF.db"):
		self.database = db
		self.chromosomes = {}
		self.log = {
			"readGFF":[],
			"readBam":[],
			"readProtein":[],
			"selectGeneModels":[],
			"writeSelected":[]
		}
		
		con = duckdb.connect(self.database)
		try: con.sql("DROP TABLE annotations CASCADE;")
		except: pass
		try: con.sql("DROP SEQUENCE seq CASCADE;")
		except: pass

		# Generate tables
		# TODO: create a view of the annotations table with de-duplicated features and their source affiliations
		self.loci_sql = """
			CREATE OR REPLACE TABLE loci (
				locus INT PRIMARY KEY DEFAULT nextval('seq'),
				start INT,
				fin INT,
				chromosome VARCHAR(50),
				strand CHAR(1) CHECK (strand IN ('+', '-')),
				selected VARCHAR
				);
			"""
		self.data_sql = '''
			CREATE OR REPLACE TABLE data (
				chromosome VARCHAR(50),
				position BIGINT
			);
			'''
		self.annotations_sql = '''
			CREATE OR REPLACE TABLE annotations (
				source VARCHAR(50),
				type VARCHAR(50),
				start INT,
				fin INT,
				gene_id VARCHAR(50),
				mrna_id VARCHAR(100),
				feature_id VARCHAR(100),
				locus INT,
				attr VARCHAR(100),
				score FLOAT,
				strand CHAR(1) CHECK (strand IN ('+', '-')),
				phase INT,
				chromosome VARCHAR(50),
				selected BOOLEAN DEFAULT FALSE
			);
			'''
		self.features_sql = """
			CREATE OR REPLACE VIEW features AS
			SELECT DISTINCT source, type, start, fin, mrna_id, locus, chromosome
			FROM annotations
			WHERE contains(type, 'CDS') OR contains(type, 'prime_UTR');
		"""
		con.sql("CREATE OR REPLACE SEQUENCE seq START 1;")
		con.sql(self.loci_sql)
		con.sql(self.data_sql)
		con.sql(self.annotations_sql)
		con.sql(self.features_sql)
		
		con.close()

	def countChromosomes(self, fasta:str):
		with open(fasta, 'r') as input:
			current_chromosome = None
			for line in input:
				line = line.strip()
				if line.startswith(">"):
					current_chromosome = line.split(" ")[0][1:]
					self.chromosomes[current_chromosome] = 0
				else:
					self.chromosomes[current_chromosome] += len(line)

	def readGFF(self, file:str, source:str):
		# Read a GFF file into the database structure.
		# Overlapping loci from distinct annotations are 
		# held at the locus level
		# TODO: GETA mrna/gene ids not properly parsed
		# TODO: could this all be ported into sql?
		con = duckdb.connect(self.database)
		
		with open(file, 'r') as gff:
			skipFlag = False
			gene_id = None
			mRNA_id = None
			locus_id = None

			num_lines = sum(1 for line in open(file,'r'))
			print(f"\nReading {file}: ", file=sys.stderr, flush=True)
			for line in tqdm(gff, total=num_lines):
				if not line:
					break
				elif line.startswith('#') | (line.strip() == ''):
					continue
				else:
					line = line.strip().split('\t')
					chr = line[0]
					typ = line[2]
					start = int(line[3])
					end = int(line[4])
					glen = end-start+1
					strand = line[6]
					attr = line[8]
					if line[7] == '.': phase = 0
					else: phase = line[7]
					if line[5] == '.': score = float(0)
					else: score = float(line[5])
				
				if typ == "gene":
					skipFlag = False
					gene_id = re.search(r"ID=([a-zA-Z0-9\.\-\_]+);{0,1}", attr)[1]
					
					# Find overlapping loci. 
					## Assumes start and end coords are from 5'-end  regardless of strand
					overlaps = con.sql(
						'SELECT locus, start, fin ' 
						'FROM loci '
						f"WHERE chromosome = '{chr}' AND strand = '{strand}' AND NOT (start > {end} OR fin < {start});"
						).df()

					if overlaps.empty:
						# Generate new loci
						con.sql(f"INSERT INTO loci (chromosome, start, fin, strand) VALUES ('{chr}',{start},{end},'{strand}');")
						locus_id = con.sql("SELECT currval('seq') AS seq;").df().loc[0,'seq']
					else:
						try:
							if max([0, min([int(overlaps['fin']), end]) - max([int(overlaps['start']), start])])/glen < 0.5: #50% overlap threshold for new loci creation
								# Generate new loci
								# troublesome spot = Chr4:10,449,966-10,456,188
								con.sql(f"INSERT INTO loci (chromosome, start, fin, strand) VALUES ('{chr}',{start},{end},'{strand}');")
								locus_id = con.sql("SELECT LAST_INSERT_ROWID();")
							else:
								# Add to existing loci
								overlaps = overlaps.iloc[0] #TODO: what if more than one overlapping locus?
								locus_id = overlaps['locus']
								
								## Update loci coordinates if necesary

								if start < overlaps['start']:
									con.sql(f"UPDATE loci SET start = {start} WHERE locus = {locus_id};")
								if end > overlaps['fin']:
									con.sql(f"UPDATE loci SET fin = {end} WHERE locus = {locus_id};")

						except:
							#TODO: how to handle loci like Chr4:10,449,966-10,456,188?
							skipFlag = True
							continue
					
					if not skipFlag:
						con.sql("INSERT INTO annotations (chromosome, source, type, start, fin, score, strand, phase, attr, gene_id, mrna_id, feature_id, locus) "
							f"VALUES ('{chr}', '{source}', '{typ}', {start}, {end}, {score}, '{strand}', {phase}, '{attr}','{gene_id}', NULL, '{gene_id}', {locus_id})")
								
				elif typ == 'ncRNA_gene':
					skipFlag = True

				elif skipFlag:
					continue

				elif typ == "mRNA":
					if "primary=False" in attr: # Skipping Mikado alternative transcripts
						skipFlag = True
						continue
					mRNA_id = re.search(r"ID=([a-zA-Z0-9\.\-\_]+);", attr)[1]
					con.sql("INSERT INTO annotations (chromosome, source, type, start, fin, score, strand, phase, attr, gene_id, mrna_id, feature_id, locus) "
							f"VALUES ('{chr}', '{source}', '{typ}', {start}, {end}, {score}, '{strand}', {phase}, '{attr}', '{gene_id}', '{mRNA_id}', '{mRNA_id}', {locus_id})")
				elif typ in ['CDS', 'exon', 'five_prime_UTR', 'three_prime_UTR']:
					feat_id = re.search(r"ID=([a-zA-Z0-9\.\-\_]+);", attr)[1]
					con.sql("INSERT INTO annotations (chromosome, source, type, start, fin, score, strand, phase, attr, gene_id, mrna_id, feature_id, locus) "
							f"VALUES ('{chr}', '{source}', '{typ}', {start}, {end}, {score}, '{strand}', {phase}, '{attr}', '{gene_id}', '{mRNA_id}', '{feat_id}', {locus_id})")
		
		# TODO: why are some whole gene models duplicated in annotations table?
		con.close()
		self.log["readGFF"].append(file)

	def readBam(self, title:str, file:str, threads=1, resume=False):
			
		wd = os.path.dirname(sys.argv[0])
		if wd == "": wd = "."

		# Compute base depth if needed
		if not os.path.exists(f"{title}.basedepth"):
			print(f"\nCalculating base depth using {threads} threads", file=sys.stderr, flush=True)
			os.system(f"{wd}/bamBaseDepth.sh {title} {file} {threads}")
		else:
			if not resume:
				print(f"\nCalculating base depth using {threads} threads", file=sys.stderr, flush=True)
				os.system(f"{wd}/bamBaseDepth.sh {title} {file} {threads}")
			else:
				print(f"\nDepth calculation skipped, {title}.basedepth exists", file=sys.stderr, flush=True)
		
		print(f"Adding {title} data to evidence.", file=sys.stderr, flush=True)
		with duckdb.connect(self.database) as con:
			con.execute("DROP INDEX IF EXISTS idx_data;")
			columns = con.execute("PRAGMA table_info('data');").df() # Check if the column exists
			if title in columns['name'].values:
				con.execute(f"ALTER TABLE data DROP COLUMN {title};") # Drop the column if it exists
			con.execute(f"ALTER TABLE data ADD COLUMN {title} FLOAT;") # Add the new column
			con.execute(f"CREATE INDEX idx_data ON data(chromosome, position);")

			# Add base depth to data table
			if con.sql("SELECT count(position) as pos from data;").fetchnumpy()["pos"][0] == 0:
				print(f"Copying from {title}.basedepth to data.", file=sys.stderr, flush=True)
				con.sql(f"COPY data FROM '{title}.basedepth' (DELIMITER ' ', HEADER TRUE, QUOTE '\"', ESCAPE '\', NULL 'NULL', AUTO_DETECT FALSE);")
			else:
				total_lines = int(subprocess.run(['wc', '-l', f"{title}.basedepth"], stdout=subprocess.PIPE).stdout.split()[0])
				chunksize = 10**8
				with tqdm(total=total_lines//chunksize, desc="Adding bam depth chunks.", unit="line") as pbar:
					for depth in pd.read_csv(f"{title}.basedepth", sep=' ', header=0, chunksize=chunksize):
						con.sql(
							"UPDATE data "
							f"SET {title} = depth.{title} "
							"FROM depth "
							"WHERE data.chromosome = depth.chromosome "
							"AND data.position = depth.position;"
							)
						pbar.update(1)

		self.log["readBam"].append(file)

	def readProtein(self, title:str, file:str, fasta:str, resume=False):
		# Read a gff formated protein alignment (such as output by miniprot) into the data table
		
		## Calculate chromosome lengths
		self.countChromosomes(fasta)

		## Write chromosomes to bed file
		if not os.path.exists("genome.bed"):
			with open("genome.bed", 'w') as output:
				for chr in self.chromosomes:
					output.write(f"{chr}\t0\t{self.chromosomes[chr]}\n")
		else:
			if not resume:
				with open("genome.bed", 'w') as output:
					for chr in self.chromosomes:
						output.write(f"{chr}\t0\t{self.chromosomes[chr]}\n")
		
		## Convert protein GFF output to bed
		if not os.path.exists(f"protein_alignments_{title}.bed"):
			os.system(f"grep CDS {file} | cut -f1,4,5 > protein_alignments_{title}.bed")
		else:
			if not resume:
				os.system(f"grep CDS {file} | cut -f1,4,5 > protein_alignments_{title}.bed")

		## Calculate depth via bedtools
		if not os.path.exists("protein_alignments_{title}.depth"):
			print(f"\nCalculating protein alignment depth.", file=sys.stderr, flush=True)
			os.system(f"bedtools coverage -d -a genome.bed -b protein_alignments_{title}.bed > protein_alignments_{title}.depth")
		else:
			if not resume:
				print(f"\nCalculating protein alignment depth.", file=sys.stderr, flush=True)
				os.system(f"bedtools coverage -d -a genome.bed -b protein_alignments_{title}.bed > protein_alignments_{title}.depth")
			else:
				print(f"\nDepth calculation skipped, protein_alignments_{title}.depth exists and {resume=}", file=sys.stderr, flush=True)

		# Add results to data table
		print("Adding protein data to evidence.", file=sys.stderr, flush=True)
		with duckdb.connect(self.database) as con:
			con.execute("DROP INDEX IF EXISTS idx_data;")
			columns = con.execute("PRAGMA table_info('data');").df() # Check if the column exists
			if title in columns['name'].values:
				con.execute(f"ALTER TABLE data DROP COLUMN {title};") # Drop the column if it exists
			con.execute(f"ALTER TABLE data ADD COLUMN {title} FLOAT;") # Add the new column
			con.execute(f"CREATE INDEX idx_data ON data(chromosome, position);")
		
		# Stream the protein_alignments.depth file in chunks (100Mbp chunks)
		total_lines = int(subprocess.run(['wc', '-l', f"protein_alignments_{title}.depth"], stdout=subprocess.PIPE).stdout.split()[0])

		chunksize = 10**8
		with tqdm(total=total_lines//chunksize, desc="Adding protein depth chunks.", unit="line") as pbar:
			for chunk in pd.read_csv(f"protein_alignments_{title}.depth", sep='\t', header=None, chunksize=chunksize):
				chunk.columns = ["chromosome", "start", "fin", "position", title]
				with duckdb.connect(self.database) as con:
					con.register("chunk", chunk)
					# Perform the SQL update
					con.sql(
						f"""
						UPDATE data 
						SET {title} = chunk.{title}
						FROM chunk
						WHERE data.chromosome = chunk.chromosome
						AND data.position = chunk.position;
						""")
				pbar.update(1)
				sys.stderr.flush()

		self.log["readProtein"].append(file)

	def selectGeneModels(self, penalty=0.01):
		#TODO: compare the annotations and select the one with the most supporting evidence
		# Idea: aggregate evidence at the feature level for each annotation 
		#       or find where the annotations diverge and see what the evidence supports?
		#       do I need to normalize the total evidence to an evidence density to not bias towards longer features?
		#       which feature is more likely to give rise to the evidence... can we use likelihoods?

		# Question: can I deal with isoforms better than just keeping all isoforms from a source?

		# Algorithm A
		# Identify the set of divergent features
		# Calculate the AED for each divergent feature set
		# Select the annotation with the lowest AED for divergent features
		
		# Algorithm B
		# Can we select individual features (regardless of source) to build a better primary transcript set
		# from the combined efforts of each annotation? 

		# Algorithm C
		# Find the transcript-wise length normalized total evidence score for each isoform
		# Select the entire annotation (source) with the highest average score among all isoforms (at each locus)
		# Log transform the evidence at each base so high evidence counts don't obscure information from low evidence counts
		# Bases with no evidence are penalized according to penalty (values between 0 and 1 are negative after log transformation)
		print("\nSelecting gene models...", file=sys.stderr, flush=True)
		con = duckdb.connect(self.database)
		try:
			con.sql("ALTER TABLE data ADD COLUMN score FLOAT;")
		except:
			con.sql("UPDATE data SET score = 0;")
		
		try:
			con.sql('DROP TABLE isoform_scores;')
			con.sql('DROP TABLE sel;')
		except:
			pass
		
		data_cols = list(con.sql("select * from data limit 1;").df().columns)
		data_cols = [col for col in data_cols if col not in ["chromosome", "position", "score"]]
		for ev in data_cols:
			# Add penalty to bases with no evidence
			con.sql('UPDATE data '
					f'SET {ev} = {penalty} '
					f'WHERE {ev} = 0;')
			con.sql('UPDATE data '
					f'SET score = score + LOG({ev}) + 1 '
				)
		
		chromosomes = [row[0] for row in con.sql('SELECT DISTINCT chromosome FROM data;').fetchall()]
		# Compute the length normalized, log transformed evidence score for each transcript

		for chrom in chromosomes:
			print(f"Processing chromosome {chrom}...", file=sys.stderr, flush=True)
			#con.sql(
			#	'CREATE TABLE isoform_scores AS '
			#	'WITH cteScores (score, mrna_id, source, locus) AS ('
			#		f'SELECT SUM(score) as score, feat.mrna_id, feat.source, feat.locus '
			#	'''FROM data, features feat
			#		WHERE data.position BETWEEN feat.start AND feat.fin
			#		GROUP BY feat.source, feat.mrna_id, feat.locus
			#	)
			#	SELECT sc.*, len.length, sc.score/len.length AS norm_score
			#	FROM cteScores sc
			#	JOIN (
			#		SELECT SUM(fin - start) as length, mrna_id
			#		FROM features
			#		GROUP BY mrna_id, source, locus) len
			#	ON sc.mrna_id = len.mrna_id
			#	''')
			print("Computing evidence scores...", file=sys.stderr)
			con.sql(f'''
				CREATE TABLE isoform_scores AS
				WITH cteScores AS (
					SELECT SUM(d.score) AS score, feat.mrna_id, feat.source, feat.locus
					FROM data d
					JOIN features feat
						ON d.position BETWEEN feat.start AND feat.fin
					WHERE d.chromosome = '{chrom}' AND feat.chromosome = '{chrom}'
					GROUP BY feat.source, feat.mrna_id, feat.locus
				)
				SELECT sc.*, len.length, sc.score / len.length AS norm_score
				FROM cteScores sc
				JOIN (
					SELECT SUM(fin - start) AS length, mrna_id
					FROM features
					WHERE chromosome = '{chrom}'
					GROUP BY mrna_id, source, locus
				) len ON sc.mrna_id = len.mrna_id;
			''')

			print("Selecting loci...", file=sys.stderr)
			# Select the annotation at each locus with maximal score 
			con.sql('''
				CREATE TABLE sel AS 
				SELECT loci.locus, loci.start, loci.fin, loci.chromosome, loci.strand, sel.selected
				FROM loci
				JOIN (
					SELECT max(sc.norm_score), sc.locus, sc.source AS selected
					FROM isoform_scores sc
					GROUP BY sc.locus, sc.source) sel
				ON loci.locus = sel.locus
				''')
			
			print("Updating loci table with selected annotations...", file=sys.stderr)
			# Update the on-disk loci table with selected annotations
			con.sql(f'''
				UPDATE loci
				SET selected = (
					SELECT sel.selected
					FROM sel
					WHERE loci.locus = sel.locus)
				WHERE loci.chromosome = '{chrom}';
				''')
			
			con.sql('DROP TABLE isoform_scores;')
			con.sql('DROP TABLE sel;')
		con.close()
		self.log["selectGeneModels"].append(penalty)
		print("Done.", file=sys.stderr, flush=True)

	def writeSelected(self, file="selected_genemodels.gff3"):
		# Parse all selected gene models into a GFF
		con =  duckdb.connect(self.database)
		con.sql("""
			COPY (
				SELECT chromosome, source, type, start, fin, score, strand, phase, attr
				FROM (
					SELECT *, MIN(start) OVER (PARTITION BY gene_id) as min_start
					FROM (
						SELECT a.chromosome, a.source, a.type, a.start, a.fin, a.score, a.strand, a.phase, a.attr, a.gene_id, row_number() OVER () AS rownum, a.locus al
						FROM annotations a
						RIGHT JOIN loci l
						ON a.locus = l.locus AND a.source = l.selected
					)
					ORDER BY min_start ASC, rownum ASC
				)
			) """
			f"TO '{file}' (FORMAT CSV, DELIMITER '\t', HEADER false)"
		)
		con.close()
		self.log["writeSelected"].append(file)
		print(f"\nSelected gene models written to: {file}", file=sys.stderr, flush=True)

	def plotEvidence(self, col:str, locus:int):
		# Plot a histogram of the locus-level evidence data for a specific column
		# TODO: fix unexpected rectangle coloring (are the transcripts actually wrong or just the coloring?)
		with duckdb.connect(self.database) as con:
			evidence = con.sql(
				f"SELECT position, {col} "
				f"FROM data d, (SELECT * FROM loci WHERE locus = {locus}) l "
				"WHERE d.position BETWEEN l.start AND l.fin;"
				).df()
			annot = con.sql(f"SELECT * from annotations WHERE locus = {locus};").df()

		positions = list(evidence.index + 1)
		feat_start = min(evidence['position'])
		maxcount = int(max(evidence[col]))
		#positions = list(map(lambda x: x - feat_start, positions))
		plt.plot(positions, evidence[col], linestyle='-')
		plt.xlabel('Genic Position')
		plt.ylabel('Evidence count')
		plt.title(f'{annot["chromosome"].unique().item()}:{min(evidence.index)}-{max(evidence.index)} => {col}')
		plt.xticks(range(0, len(positions), 500), )
		plt.yticks(range(int(-1*maxcount/10), int(maxcount), int(maxcount/10))) # TODO: ensure third argument does not round to 0

		isoforms = [i for i in list(annot['mrna_id'].unique()) if i]
		num_annotations = len(isoforms)
		anot_height = -1*(maxcount/(6*num_annotations))

		rectangles = []
		row = 0
		for isoform in isoforms:
			for _,feat in annot[annot['mrna_id'] == isoform].iterrows():
				if (feat['type'] == 'five_prime_UTR') | (feat['type'] == 'three_prime_UTR'):
					feat_color = 'green'
				elif feat['type'] == 'CDS':
					feat_color = 'lightblue'
				#elif feat['type'] == 'exon': # this may overwrite 3'UTR in some cases
				#	feat_color = 'blue'
				else:
					feat_color = None
				if feat_color:
					plt.text(0, 
						(row+(anot_height/2)), 
						feat['source'], 
						fontsize = 8,
						verticalalignment='center',
						clip_on = True) 
					rectangles.append(
						Rectangle((feat['start']-feat_start, row),
									feat['fin']-feat['start'],
									anot_height,
									color=feat_color, 
									alpha=1,
									ec='black'))
			row += anot_height
			
		for rectangle in rectangles:
			plt.gca().add_patch(rectangle)

		custom_lines = [Line2D([0], [0], color='green', lw=5),
						Line2D([0], [0], color='lightblue', lw=5)]
		plt.legend(custom_lines, ['UTR', 'CDS'])
			
		for rectangle in rectangles:
			plt.gca().add_patch(rectangle)

		plt.grid(True)
		plt.show()

	def runCombine(self, gffs:list, bams:list, prots:list, fasta:str, output:str, pkl="EvidentialCombine.pkl", resume=False, sw=False):
		if not sw:
			for gff in gffs:
				file, title = gff.split(":")
				if resume & (file in self.log["readGFF"]): 
					print(f"Skipping {gff} because {resume=}", file=sys.stderr, flush=True)
					continue
				self.readGFF(file, title)
				with open(pkl, "wb") as p: 
					pickle.dump(self,p)

			for bam in bams:
				file, title = bam.split(":")
				if resume & (file in self.log["readBam"]):
					print(f"Skipping {bam} because {resume=}", file=sys.stderr, flush=True)
					continue
				self.readBam(title, file, resume=resume)
				with open(pkl, "wb") as p: 
					pickle.dump(self,p)
			
			for prot in prots:
				file, title = prot.split(":")
				if resume & (file in self.log["readProtein"]):
					print(f"Skipping {prot} because {resume=}", file=sys.stderr, flush=True)
					continue
				self.readProtein(title, file, fasta, resume)
				with open(pkl, "wb") as p: 
					pickle.dump(self,p)
		
		self.selectGeneModels()
		with open(pkl, "wb") as p: 
				pickle.dump(self,p)
		
		self.writeSelected(output)
		with open(pkl, "wb") as p: 
				pickle.dump(self,p)

	# TODO: Write annotation comparison method (implement annotation edit distance)

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument("-g", "--gff_files", default=[], nargs="+", help="GFF files to be combined. Provided titles following colon (gff:title)")
	parser.add_argument("-b", "--bam_files", default=[], nargs="+", help="RNAseq alignent files (BAM). Provided titles following colon (bam:title)")
	parser.add_argument("-p", "--prot_files", default=[], nargs="+", help="Protein alignment files (gff). Provided titles following colon (protein:title)")
	parser.add_argument("-f", "--genome_fasta", type=str, help="Genome fasta file. Required for protein alignment reading.")
	parser.add_argument("-i", "--existing", type=str, default="", help="Pickle file of existing run")
	parser.add_argument("-o", "--output", type=str, default="selected_gene_models.gff3", help="Name of output file gff")
	parser.add_argument("--pkl", type=str, default="EvidentialCombine.pkl", help="Name of output pickled genome annotation obect")
	parser.add_argument("-d", "--database", type=str, default="CombineGFF.db", help="Name of DuckDB database to use")
	parser.add_argument("--resume", action="store_true", help="Resume from last successful point, otherwise start from beginning")
	parser.add_argument("--selectandwrite", action="store_true", help="Only select and write gene models")
	args = parser.parse_args()

	if args.existing:
		with open(args.existing, 'rb') as p:
			ga = pickle.load(p)
		print(f"Running from existing: {args.existing}...database={ga.database}, resume={args.resume}")
	else:
		print(f"Running from scratch...creating new {args.existing} with database: {args.database}")
		ga = genome(db=args.database)
	
	ga.runCombine(
		args.gff_files, 
		args.bam_files, 
		args.prot_files, 
		args.genome_fasta, 
		args.output, 
		pkl=args.pkl, 
		resume=args.resume, 
		sw=args.selectandwrite
	)
	
	print("Complete.")