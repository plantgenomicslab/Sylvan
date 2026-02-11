# Sylvan Genome Annotation

Sylvan is a comprehensive genome annotation pipeline that combines EVM/PASA, GETA, and Helixer with semi-supervised random forest filtering for generating high-quality gene models from raw genome assemblies.

![Sylvan Workflow](https://github.com/plantgenomicslab/Sylvan/blob/main/docs/images/sylvan_workflow.jpg?raw=true)

## Pipeline Design

The Sylvan pipeline consists of two main phases — annotation and filtration — with configurable modules that process evidence from multiple sources and combine them into a unified gene model. The following describes the available tools and modules. Users configure which components to enable and how to parameterize them via `config_annotate.yml` and `config_filter.yml`.

![Pipeline DAG](https://github.com/plantgenomicslab/Sylvan/blob/main/docs/images/rulegraph.png?raw=true)

---

### Phase 1: Annotate

The annotation phase generates gene models by integrating multiple configurable evidence sources.

#### Evidence Generation

- **Repeat Masking**
  - Runs RepeatMasker with a user-specified species library (e.g. `Embryophyta`, `Viridiplantae`, `Metazoa` — configured via `geta.RM_species`)
  - Can optionally run RepeatModeler for de novo repeat identification
  - Supports user-supplied custom repeat libraries (e.g. from EDTA, configured via `geta.RM_lib`)

- **RNA-seq Processing**
  - Quality-trims reads with fastp
  - Aligns reads via STAR (default) or HiSat2 (alternative pathway — both are available in the pipeline; the active pathway depends on the Snakemake rule graph)
  - Assembles transcripts with StringTie and PsiCLASS
  - Optionally performs de novo transcript assembly with SPAdes + Evigene clustering
  - Refines and clusters transcripts with PASA

- **Protein Homology** (sequential pipeline)
  - Miniprot performs fast protein-to-genome alignment to identify candidate gene regions
  - GeneWise refines gene structures on Miniprot-identified regions
  - GMAP provides exonerate-style exon-level alignments

- **Ab Initio Prediction**
  - Helixer: deep learning–based gene prediction (optionally GPU-accelerated; model selected via `helixer_model` — `land_plant`, `vertebrate`, `invertebrate`, or `fungi`)
  - Augustus: HMM-based prediction, either trained de novo on the target genome or initialized from an existing species model (via `augustus_start_from`), or skipped entirely if a pre-trained model is supplied (via `use_augustus`)

- **Liftover**
  - LiftOff transfers annotations from one or more neighbor species (configured via `liftoff.neighbor_gff` and `liftoff.neighbor_fasta`)

#### Evidence Combination

- **GETA Pipeline**
  - TransDecoder predicts ORFs from assembled transcripts
  - Gene models are combined and filtered; repeat-overlapping genes are removed

- **Portcullis**
  - Filters splice junctions from transcript evidence

- **EvidenceModeler (EVM)**
  - Integrates all evidence sources using configurable weights (`evm_weights.txt`)
  - Generates consensus gene models
  - Genome is partitioned into overlapping segments for parallel execution (partition count configured via `num_evm_files`)

- **PASA Post-processing**
  - PASA operates at two stages in the pipeline: (1) initial transcript assembly and clustering before EVM, and (2) post-EVM refinement for UTR addition and alternative isoform incorporation

- **AGAT**
  - Final GFF3 format cleaning and validation

**Output:** `results/complete_draft.gff3`

---

### Phase 2: Semi-Supervised Random Forest Filter

The filter phase computes additional evidence features for each gene model and applies a semi-supervised random forest classifier to separate high-quality genes from spurious predictions.

#### Feature Generation

The following features are computed for every gene model in the draft annotation:

- **PfamScan** — identifies conserved protein domains using the Pfam-A HMM database
- **RSEM** — quantifies transcript expression (TPM) from re-aligned RNA-seq reads; bedtools computes read coverage
- **BLASTp (homolog)** — measures similarity to a user-supplied protein database (parallelized across 20 split peptide files)
- **BLASTp (RexDB)** — measures similarity to a repeat element protein database (e.g. RepeatExplorer Viridiplantae)
- **Ab initio overlap** — computes the fraction of each gene model overlapping with Augustus predictions, Helixer predictions, and RepeatMasker annotations
- **lncDC** — classifies transcripts as protein-coding or long non-coding RNA using an XGBoost model with plant-specific pre-trained parameters
- **BUSCO** — identifies conserved single-copy orthologs (*used only to monitor the filtration process; not used as a classifier feature*)

#### Semi-supervised Classification

1. **Initial gene set selection**: A data-driven heuristic selects high-confidence positive genes (strong homolog/Pfam/expression evidence) and high-confidence negative genes (repeat-like, no expression) using configurable cutoff thresholds (TPM, coverage, BLAST identity/coverage, repeat overlap)
2. **Random forest training**: A binary classifier is trained on the initial gene set
3. **Iterative refinement**: High-confidence predictions (above the `--recycle` threshold, default 0.95) are added back to the training set, and the model is retrained. This repeats for up to `--max-iter` iterations (default 5) or until convergence

**Output files:**
- `results/FILTER/filter.gff3` — Kept gene models
- `results/FILTER/discard.gff3` — Discarded gene models
- `results/FILTER/data.tsv` — Feature matrix used by random forest
- `results/FILTER/keep_data.tsv` — Evidence data for kept genes
- `results/FILTER/discard_data.tsv` — Evidence data for discarded genes
- `results/FILTER/{prefix}.cdna` — Extracted transcript sequences
- `results/FILTER/{prefix}.pep` — Extracted peptide sequences

---

### Alternative: Score-based Filter

An alternative scoring pipeline (`Snakefile_filter_score`) uses logistic regression and random forest scoring with pseudo-labels instead of the iterative semi-supervised approach. This requires the same feature generation outputs and produces:
- `results/FILTER/scores.csv` — Per-gene scores and features
- `results/FILTER/scores.metrics.txt` — AUC/PR/F1 and chosen thresholds

```bash
export SYLVAN_FILTER_CONFIG="toydata/config/config_filter.yml"
./bin/filter_score_toydata.sh
```

---

## Features

- **Multi-evidence integration**: RNA-seq, protein homology, neighbor species annotations
- **Dual RNA-seq alignment pathways**: STAR and HiSat2 with StringTie/PsiCLASS
- **Multiple ab initio predictors**: Helixer (GPU-accelerated), Augustus
- **Semi-supervised filtering**: Random forest-based spurious gene removal
- **Score-based filtering**: Alternative logistic regression + random forest scoring pipeline
- **HPC-ready**: SLURM cluster support with Singularity containers
- **TidyGFF**: Format annotations for public distribution
- **Cleanup utility**: Remove intermediate files after pipeline completion

## Quick Start

```bash
# 1. Install environment
conda create -n sylvan -c conda-forge -c bioconda python=3.11 snakemake=7 -y
conda activate sylvan

# 2. Download Singularity image
singularity pull --arch amd64 sylvan.sif library://wyim/sylvan/sylvan:latest

# 3. Clone repository
git clone https://github.com/plantgenomicslab/Sylvan.git
cd Sylvan

# 4. Run with toy data (dry-run first)
snakemake -n --snakefile bin/Snakefile_annotate
./bin/annotate_toydata.sh
```

The toy data experiment uses *A. thaliana* chromosome 4 with 12 paired-end RNA-seq samples, 3 neighbor species, and the `land_plant` Helixer model. For a detailed walkthrough of this experiment, see the **[Wiki](Wiki.md)**.

Helper script:
- `bin/generate_cluster_from_config.py`: optionally regenerate per-rule SLURM defaults within `config_annotate.yml` — keeps resource requests in sync with the pipeline's threads/memory.

## Installation

### Requirements

- Linux (tested on CentOS/RHEL)
- Singularity 3.x+
- Conda/Mamba
- SLURM (for cluster execution)
- Git LFS (for toy data)

### Dependencies

Most bioinformatics tools (STAR, Augustus, GeneWise, PASA, EVM, BLAST, BUSCO, etc.) are bundled inside the Singularity container. The host environment needs:

| Package | Purpose |
|---------|---------|
| Python 3.10+ | Pipeline orchestration |
| Snakemake 7 | Workflow engine |
| pandas | Data manipulation |
| scikit-learn | Random forest classifier |
| NumPy | Numerical operations |
| PyYAML | Config parsing |
| rich | Logging (optional) |

Perl and R scripts (`fillingEndsOfGeneModels.pl`, `filter_distributions.R`) run inside the Singularity container and do not require host installation.

### Setup

```bash
# Create conda environment
conda create -n sylvan -c conda-forge -c bioconda python=3.11 snakemake=7 -y
conda activate sylvan

# Download Singularity image
singularity pull --arch amd64 sylvan.sif library://wyim/sylvan/sylvan:latest

# Clone repository (with Git LFS for toy data)
git lfs install
git clone https://github.com/plantgenomicslab/Sylvan.git
```

### Build from source (optional)

```bash
cd Sylvan/singularity
sudo singularity build sylvan.sif Sylvan.def
```

## Running the Annotate Phase

This section describes the inputs, configuration, and commands needed to run the annotation pipeline on your data.

### Input Requirements

| Input | Description | Config Field |
|-------|-------------|--------------|
| Genome assembly | FASTA file (`.fa`, `.fasta`, `.fa.gz`, `.fasta.gz`) | `genome` |
| RNA-seq data | Paired-end gzipped FASTQ files (`*_1.fastq.gz`/`*_2.fastq.gz` or `*_R1.fastq.gz`/`*_R2.fastq.gz`) in a folder | `rna_seq` |
| Protein sequences | FASTA from UniProt, OrthoDB, etc. (comma-separated for multiple files) | `proteins` |
| Neighbor species | Directories containing GFF3 and genome FASTA files (one file per species) | `liftoff.neighbor_gff`, `liftoff.neighbor_fasta` |
| Repeat library | EDTA output (`.TElib.fa`) | `geta.RM_lib` |
| Singularity image | Path to `sylvan.sif` | `singularity` |

### Running the Pipeline

```bash
# Set config (required)
export SYLVAN_CONFIG="toydata/config/config_annotate.yml"

# Dry run
snakemake -n --snakefile bin/Snakefile_annotate

# Submit to SLURM
sbatch -A [account] -p [partition] -c 1 --mem=1g \
  -J annotate -o annotate.out -e annotate.err \
  --wrap="./bin/annotate_toydata.sh"

# Or run directly
./bin/annotate_toydata.sh

```

**Output:** `results/complete_draft.gff3`

---

## Running the Filter Phase

This section describes the inputs and commands for the filter pipeline. All inputs below are specified in `config_filter.yml`.

### Input Requirements

| Input | Description | Config Field |
|-------|-------------|--------------|
| Annotated GFF | Output from Annotate phase (`results/complete_draft.gff3`) | `anot_gff` |
| Genome | Same as Annotate phase | `genome` |
| RNA-seq data | Same as Annotate phase | `rna_seq` |
| Protein sequences | Same as Annotate phase | `protein` |
| Augustus GFF | Augustus predictions (`results/GETA/Augustus/augustus.gff3`) | `augustus_gff` |
| Helixer GFF | Helixer predictions (`results/AB_INITIO/Helixer/helixer.gff3`) | `helixer_gff` |
| Repeat GFF | RepeatMasker output (`results/GETA/RepeatMasker/genome.repeat.gff3`) | `repeat_gff` |
| HmmDB | Pfam database directory (default: `/usr/local/src` inside container) | `HmmDB` |
| RexDB | RepeatExplorer protein DB (e.g. `Viridiplantae_v4.0.fasta` from [rexdb](https://github.com/repeatexplorer/rexdb)) | `RexDB` |
| BUSCO lineage | e.g., `eudicots_odb10` | `busco_lin` |
| Chromosome regex | Regex to match chromosome prefixes (e.g. `(^Chr)\|(^chr)\|(^LG)`) | `chrom_regex` |

**Filter cutoff thresholds** (in `config_filter.yml` under `Cutoff`):

| Parameter | Description | Default |
|-----------|-------------|---------|
| `tpm` | TPM threshold for initial gene selection | 3 |
| `rsem_cov` | RNA-seq coverage threshold | 0.5 |
| `blast_pident` | BLASTp percent identity threshold | 0.6 |
| `blast_qcovs` | BLASTp query coverage threshold | 0.6 |
| `rex_pident` | RexDB percent identity threshold | 0.6 |
| `rex_qcovs` | RexDB query coverage threshold | 0.6 |
| `helixer_cov` | Helixer overlap coverage threshold | 0.8 |
| `augustus_cov` | Augustus overlap coverage threshold | 0.8 |
| `repeat_cov` | Repeat overlap coverage threshold | 0.5 |

### Running the Pipeline

```bash
# Set config (required)
export SYLVAN_FILTER_CONFIG="toydata/config/config_filter.yml"

# Dry run
snakemake -n --snakefile bin/Snakefile_filter

# Submit to SLURM
sbatch -A [account] -p [partition] -c 1 --mem=4g \
  -J filter -o filter.out -e filter.err \
  --wrap="./bin/filter_toydata.sh"

# Or run directly
./bin/filter_toydata.sh
```

**Output:** `results/FILTER/filter.gff3`

### Feature Importance Test

Reviewers often ask for an ablation study of the semi-supervised filter. After a
filter run completes (which produces `FILTER/data.tsv`), launch the automated
leave-one-feature-out test:

```bash
python bin/filter_feature_importance.py FILTER/data.tsv results/busco/full_table.tsv \
  --output-table FILTER/feature_importance.tsv
```

The script reuses `Filter.semiSupRandomForest`, trains a baseline model with all
features, and then retrains while removing each feature individually. The final
out-of-bag error deltas are written to `FILTER/feature_importance.tsv` (and
`FILTER/feature_importance.json`). Use `--features` to restrict the analysis to a
subset of columns or `--ignore` to drop metadata columns that should never be
used as predictors.

## Configuration

Sylvan uses several configuration files:

| File | Purpose |
|------|---------|
| `config_annotate.yml` | **Pipeline options and SLURM resources**: input paths, species parameters, tool settings, plus per-rule CPU/memory/partition allocation |
| `config_filter.yml` | **Filter options**: input paths, cutoff thresholds, thread allocation |
| `evm_weights.txt` | **EVM evidence weights**: priority of each evidence source |
| `config/plant.yaml` | **Mikado scoring**: transcript selection parameters (plant-specific defaults provided) |

### Pipeline Config (`config_annotate.yml`)

Contains:
- Input file paths (genome, RNA-seq, proteins, neighbor species)
- Species-specific settings (Helixer model, Augustus species)
- Tool parameters (max intron length, EVM weights)
- Output prefix and directories
- SLURM resource allocation (`__default__` section with account, partition, memory, ncpus, time, plus per-rule overrides)

### EVM Weights (`evm_weights.txt`)

Controls how EvidenceModeler prioritizes different evidence sources. Higher weights give more influence. Example (from toy data):

```
ABINITIO_PREDICTION  AUGUSTUS     7
ABINITIO_PREDICTION  Helixer     3
OTHER_PREDICTION     Liftoff     2
OTHER_PREDICTION     GETA        5
OTHER_PREDICTION     Genewise    2
TRANSCRIPT           assembler-pasa.sqlite  10
TRANSCRIPT           StringTie   1
TRANSCRIPT           PsiClass    1
PROTEIN              GeneWise    2
PROTEIN              miniprot    2
```

Adjust weights based on the quality of each evidence type for your organism. PASA transcripts (weight 10) typically have the highest weight as they represent direct transcript evidence.

### Environment Variables

| Variable | Phase | Description |
|----------|-------|-------------|
| `SYLVAN_CONFIG` | Annotate | Path to `config_annotate.yml` (default: `config_annotate.yml` in cwd) |
| `SYLVAN_FILTER_CONFIG` | Filter | Path to `config_filter.yml` (default: `config_filter.yml` in cwd) |
| `SYLVAN_RESULTS_DIR` | Annotate | Override results output directory (default: `$(pwd)/results/`) |
| `TMPDIR` | Both | Temporary directory — **critical on HPC** (see below) |
| `SLURM_TMPDIR` | Both | Should match `TMPDIR` |
| `SINGULARITY_BIND` | Both | Bind additional host paths into container |

**Why `TMPDIR` matters:** Many HPC nodes mount `/tmp` as `tmpfs` (RAM-backed). Large temporary files from STAR, RepeatMasker, or Augustus can exhaust memory, causing cryptic segmentation faults or "no space left on device" errors. Always set `TMPDIR` to disk-backed project storage:

```bash
mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"
export SLURM_TMPDIR="$TMPDIR"
```

### Using Custom Config Location

```bash
# For toydata
export SYLVAN_CONFIG="toydata/config/config_annotate.yml"

# For custom project
export SYLVAN_CONFIG="/path/to/my_config.yml"
```

This is required for any Snakemake command (dry-run, unlock, etc.):

```bash
export SYLVAN_CONFIG="toydata/config/config_annotate.yml"
snakemake -n --snakefile bin/Snakefile_annotate  # dry-run
snakemake --unlock --snakefile bin/Snakefile_annotate  # unlock
```

### Key Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `prefix` | Output file prefix | `my_species` |
| `helixer_model` | `land_plant`, `vertebrate`, `invertebrate`, `fungi` | `land_plant` |
| `helixer_subseq` | 64152 (plants), 21384 (fungi), 213840 (vertebrates) | `64152` |
| `augustus_species` | Augustus species name for training | `arabidopsis` |
| `augustus_start_from` | Start Augustus training from an existing species model (skips de novo training if close match available) | `arabidopsis` |
| `use_augustus` | Use a pre-trained Augustus species without re-training (set to species name, or `placeholder` to train fresh) | `placeholder` |
| `num_evm_files` | Number of parallel EVM partitions (more = faster but more SLURM jobs) | `126` |
| `geta.RM_species` | RepeatMasker species database (e.g. `Embryophyta`, `Viridiplantae`, `Metazoa`) | `Embryophyta` |

### Helixer GPU Configuration

Helixer benefits significantly from GPU acceleration (~10x speedup). To use a separate GPU partition, add the following per-rule override in `config_annotate.yml`:

```yaml
helixer:
  ncpus: 4
  memory: 32g
  account: your-gpu-account      # GPU-specific billing account
  partition: your-gpu-partition   # GPU partition name
```

### SLURM Configuration

Find your SLURM account and partition:
```bash
# Show your accounts and partitions
sacctmgr show user "$USER" withassoc format=Account,Partition -nP

# List all available partitions
sinfo -s

# Show partition details (time limits, nodes, etc.)
sinfo -o "%P %l %D %c %m"
```

Set in `config_annotate.yml` under the `__default__` section:
```yaml
__default__:
  account: your-account
  partition: your-partition
  time: "3-00:00:00"
  memory: 8g
  ncpus: 1
```

## Useful Commands

```bash
# Force rerun all
./bin/annotate.sh --forceall

# Rerun specific rule
./bin/annotate.sh --forcerun helixer

# Rerun incomplete jobs (jobs that started but didn't finish)
./bin/rerun-incomplete.sh

# Generate report after completion
snakemake --report report.html --snakefile bin/Snakefile_annotate

# Unlock after interruption
./bin/annotate.sh --unlock

# Clean up intermediate files (run after BOTH phases complete)
./bin/cleanup.sh
```

### Cleanup

`bin/cleanup.sh` removes intermediate files generated during the annotation phase while preserving:
- Final outputs (`complete_draft.gff3`, `filter.gff3`)
- Log files (`results/logs/`)
- Configuration files
- Filter phase outputs (`FILTER/`)

Run this only after both annotation and filter phases have completed successfully.

## Output

All outputs are organized under `results/`:

```
results/
├── complete_draft.gff3          # Annotate phase final output
│
├── AB_INITIO/
│   └── Helixer/                 # Helixer predictions
│
├── GETA/
│   ├── RepeatMasker/            # Repeat masking results
│   ├── Augustus/                # Augustus predictions
│   ├── transcript/              # TransDecoder results
│   ├── homolog/                 # Protein alignments (Miniprot → GeneWise)
│   └── CombineGeneModels/       # GETA gene models
│
├── LIFTOVER/
│   └── LiftOff/                 # Neighbor species liftover
│
├── TRANSCRIPT/
│   ├── PASA/                    # PASA assemblies
│   ├── spades/                  # De novo assembly
│   └── evigene/                 # Evigene transcript clustering
│
├── PROTEIN/                     # Merged protein alignments
│
├── EVM/                         # EvidenceModeler output
│
├── FILTER/
│   ├── filter.gff3              # Kept gene models
│   ├── discard.gff3             # Discarded gene models
│   ├── data.tsv                 # Feature matrix (input to random forest)
│   ├── keep_data.tsv            # Evidence data for kept genes
│   ├── discard_data.tsv         # Evidence data for discarded genes
│   ├── {prefix}.cdna            # Transcript sequences
│   ├── {prefix}.pep             # Peptide sequences
│   ├── STAR/                    # RNA-seq realignment for filter
│   ├── rsem_outdir/             # RSEM quantification
│   ├── splitPep/                # Parallelized BLAST inputs
│   ├── busco_*/                 # BUSCO results (monitoring only)
│   └── lncrna_predict.csv       # lncDC predictions
│
├── config/                      # Runtime config copies
│
└── logs/                        # SLURM job logs
```

## Formatting Output

Use TidyGFF to prepare annotations for public distribution:

```bash
singularity exec sylvan.sif python bin/TidyGFF.py \
  MySpecies results/FILTER/filter.gff3 \
  --out MySpecies_v1.0 --splice-name t --justify 5 --sort \
  --chrom-regex "^Chr" --source Sylvan
```

**TidyGFF options:**

| Option | Description |
|--------|-------------|
| `pre` (positional) | Prefix for gene IDs (e.g. `Ath` produces `Ath01G000010`) |
| `gff` (positional) | Input GFF3 file |
| `--out` | Output file basename (produces `.gff3`, `.cdna`, `.pep` files) |
| `--splice-name` | Splice variant label style (e.g. `t` → `mRNA1.t1`, `mRNA1.t2`) |
| `--justify` | Number of digits in gene IDs (default: 8) |
| `--sort` | Sort output by chromosome and start coordinate |
| `--chrom-regex` | Regex for chromosome prefixes (auto-detects `Chr`, `chr`, `LG`, `Ch`, `^\d`) |
| `--contig-regex` | Regex for contig/scaffold naming (e.g. `HiC_scaffold_(\d+$),Scaf`) |
| `--source` | Value for GFF column 2 (e.g. `Sylvan`) |
| `--remove-names` | Remove Name attributes from GFF |
| `--no-chrom-id` | Do not number gene IDs by chromosome |

## Troubleshooting

### Check logs

```bash
# Find recent errors
ls -lt results/logs/*.err | head -10
grep -l 'Error\|Traceback' results/logs/*.err

# View specific log
cat results/logs/{rule}_{wildcards}.err
```

### Common Issues

| Issue | Solution |
|-------|----------|
| Out of memory | Increase `memory` in `config_annotate.yml` for the rule |
| `No space left on device` | `TMPDIR` is on tmpfs or quota exceeded — set `TMPDIR` to project storage |
| `Segmentation fault` | Often caused by tmpfs exhaustion — set `TMPDIR` to disk-backed storage |
| File not found (Singularity) | Path not bound in container — add to `SINGULARITY_BIND` |
| `Permission denied` in container | Check directory permissions, ensure path is bound |
| SLURM account error | Use `account` (billing account), not username |
| LFS files not downloaded | Run `git lfs pull`; verify with `ls -la toydata/` (files should be > 200 bytes) |
| Augustus training fails | Needs minimum ~500 training genes; use `augustus_start_from` with a close species |
| Job timeout | Increase `time` in `config_annotate.yml` for the rule |
| Variables not in SLURM job | Add `#SBATCH --export=ALL` or explicitly export in submit script |

### Memory Guidelines

- General recommendation: **4GB per thread**
- Example: 48 threads = 192g memory
- `ncpus` and `threads` should match in `config_annotate.yml`
- Some rules need more: `mergeSTAR` may require ~18GB per thread for large datasets
- Check `df -h $TMPDIR` to ensure temp storage is on real disk, not tmpfs

## No HPC?

Deploy a SLURM cluster on Google Cloud: [Cloud Cluster Toolkit](https://docs.cloud.google.com/cluster-toolkit/docs/quickstarts/slurm-cluster)

## Citation

> Sylvan: A comprehensive genome annotation pipeline. *Under review.*

## License

MIT License - see [LICENSE](LICENSE)

## Contact

Issues: https://github.com/plantgenomicslab/Sylvan/issues
