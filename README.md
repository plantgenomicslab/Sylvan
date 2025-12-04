# Sylvan Genome Annotation

Sylvan is a comprehensive genome annotation pipeline that combines EVM/PASA, GETA, and Helixer with semi-supervised random forest filtering for generating high-quality gene models from raw genome assemblies.

![Sylvan Workflow](https://github.com/plantgenomicslab/Sylvan/blob/main/docs/images/sylvan_workflow.jpg?raw=true)

## Pipeline Workflow

The Sylvan pipeline consists of two main phases, annotation and filtration, with interconnected modules that process evidence from multiple sources and combine them into a unified gene model.

![Pipeline DAG](https://github.com/plantgenomicslab/Sylvan/blob/main/docs/images/rulegraph.png?raw=true)

---

## Phase 1: Annotate

The annotation phase generates gene models from multiple evidence sources.

### Evidence Generation
- **Repeat Masking**
  - RepeatMasker with species-specific libraries
  - RepeatModeler for de novo repeat identification
  - Custom repeat library support (EDTA)

- **RNA-seq Processing**
  - Quality control with fastp
  - Alignment with STAR
  - Transcript assembly with StringTie and PsiCLASS
  - PASA refinement and clustering

- **Protein Homology**
  - Miniprot for fast protein-to-genome alignment
  - GeneWise for refined gene structure prediction
  - GMAP for exonerate-style alignments

- **Ab Initio Prediction**
  - Helixer (deep learning-based)
  - Augustus (HMM-based, trained on your data)

- **Liftover**
  - LiftOff for transferring annotations from neighbor species

### Evidence Combination
- **GETA Pipeline**
  - TransDecoder for ORF prediction
  - Gene model combination and filtering

- **EvidenceModeler (EVM)**
  - Weighted evidence integration
  - Consensus gene model generation

- **PASA Update**
  - UTR addition and refinement
  - Alternative isoform incorporation

**Output:** `results/complete_draft.gff3`

---

## Phase 2: Semi-Supervised Random Forest Filter

The filter phase refines and validates the annotation using additional evidence.

### Data generation
- **PfamScan**
  - Identification of conserved protein domains
- **RSEM**
  - Transcript quantification
  - RNAseq data coverage
- **BLAST**
  - Similarity to protein database
  - Similarity to repeat element database
- **lncDC**
  - Classification of long non-coding RNAs
- **BUSCO**
  - Identify conserved gene models
  - *Used only to monitor filtration process*

### Semi-supervised classification
- **Select high confidence genes**
  - Data-driven heuristic to select intial gene set
  - Select both true and spurious genes

- **Classification**
  - Train random forest binary classifier on intial gene set
  - Iteratively update gene set from predictions and re-train

**Output:** `results/FILTER/filter.gff3`

---

## Features

- **Multi-evidence integration**: RNA-seq, protein homology, neighbor species annotations
- **Multiple ab initio predictors**: Helixer, Augustus
- **Semi-supervised filtering**: Random forest-based spurious gene removal
- **HPC-ready**: SLURM cluster support with Singularity containers
- **TidyGFF**: Format annotations for public distribution

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

Helper script:
- `bin/generate_cluster_from_config.py`: regenerate `cluster_annotate.yml` from `config_annotate.yml` - used to keep SLURM resource requests in sync with the pipeline's threads/memory.

For a detailed tutorial with toy data, see the **[Wiki](Wiki.md)**.

## Installation

### Requirements

- Linux (tested on CentOS/RHEL)
- Singularity 3.x+
- Conda/Mamba
- SLURM (for cluster execution)
- Git LFS (for toy data)

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

## Annotate Phase

### Input Requirements

| Input | Description | Config Field |
|-------|-------------|--------------|
| Genome assembly | FASTA file (`.fa`, `.fasta`, `.fa.gz`, `.fasta.gz`) | `genome` |
| RNA-seq data | Gzipped FASTQ files in a folder | `rna_seq` |
| Protein sequences | FASTA from UniProt, OrthoDB, etc. | `proteins` |
| Neighbor species | GFF3 + genome FASTA files | `liftoff.neighbor_gff`, `liftoff.neighbor_fasta` |
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
  --wrap="./bin/annotate_toydata.sh
"

# Or run directly
./bin/annotate_toydata.sh

```

**Output:** `results/complete_draft.gff3`

---

## Filter Phase

### Input Requirements

| Input | Description | Config Field |
|-------|-------------|--------------|
| Annotated GFF | Output from Annotate phase | `anot_gff` |
| Genome | Same as Annotate phase | `genome` |
| RNA-seq data | Same as Annotate phase | `rna_seq` |
| Protein sequences | Same as Annotate phase | `protein` |
| Augustus GFF | Augustus predictions | `augustus_gff` |
| Helixer GFF | Helixer predictions | `helixer_gff` |
| Repeat GFF | RepeatMasker output | `repeat_gff` |
| HmmDB | Pfam database directory | `HmmDB` |
| RexDB | Plant repeat database (e.g. Viridiplantae_v4.0.fasta) | `RexDB` |
| BUSCO lineage | e.g., `eudicots_odb10` | `busco_lin` |

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

## Configuration

Sylvan uses two separate configuration files:

| File | Purpose |
|------|---------|
| `config_annotate.yml` | **Pipeline options**: input paths, species parameters, tool settings |
| `cluster_annotate.yml` | **SLURM resources**: CPU, memory, partition for each rule |

### Pipeline Config (`config_annotate.yml`)

Contains:
- Input file paths (genome, RNA-seq, proteins, neighbor species)
- Species-specific settings (Helixer model, Augustus species)
- Tool parameters (max intron length, EVM weights)
- Output prefix and directories

### Cluster Config (`cluster_annotate.yml`)

Contains SLURM resource allocation organized by pipeline phase:

```yaml
################################################################################
#                           ANNOTATE PHASE
################################################################################

#===============================================================================
# Genome Preparation
#===============================================================================
prepareGenome:
  ncpus: 1
  memory: 4g

#===============================================================================
# Repeat Masking (GETA)
#===============================================================================
RepeatMasker_species:
  ncpus: 4
  threads: 4
  memory: 16g
# ... more rules ...

#===============================================================================
# EVM - Evidence Modeler
#===============================================================================
runEVM:
  ncpus: 1
  memory: 8g

################################################################################
#                           FILTER PHASE
################################################################################

#===============================================================================
# Mikado - Transcript Selection
#===============================================================================
mikadoPick:
  ncpus: 4
  memory: 16g
```

This separation allows you to reuse the same pipeline config across different clusters by only changing the cluster config.

### Using Custom Config Location

Set the `SYLVAN_CONFIG` environment variable to use a config file in a different location:

```bash
# For toydata
export SYLVAN_CONFIG="toydata/config/config_annotate.yml"

# For custom project
export SYLVAN_CONFIG="/path/to/my_config.yml"

# The cluster config is auto-derived (cluster_annotate.yml in same directory)
# Or set explicitly:
export SYLVAN_CLUSTER_CONFIG="/path/to/my_cluster.yml"
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
| `augustus_species` | Augustus species or custom name | `arabidopsis` |

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

Set in `cluster_annotate.yml`:
```yaml
__default__:
  account: your-account
  partition: your-partition
```

## Useful Commands

```bash
# Force rerun all
./bin/annotate.sh --forceall

# Rerun specific rule
./bin/annotate.sh --forcerun helixer

# Generate report after completion
snakemake --report report.html --snakefile bin/Snakefile_annotate

# Unlock after interruption
./bin/annotate.sh --unlock
```

## Output

All outputs are organized under `results/`:

```
results/
├── complete_draft.gff3          # Annotate phase output
│
├── AB_INITIO/
│   └── Helixer/                 # Helixer predictions
│
├── GETA/
│   ├── RepeatMasker/            # Repeat masking results
│   ├── Augustus/                # Augustus predictions
│   ├── transcript/              # TransDecoder results
│   ├── homolog/                 # Protein alignments
│   └── CombineGeneModels/       # GETA gene models
│
├── LIFTOVER/
│   └── LiftOff/                 # Neighbor species liftover
│
├── TRANSCRIPT/
│   ├── PASA/                    # PASA assemblies
│   ├── spades/                  # De novo assembly
│   └── evigene/                 # Evigene clustering
│
├── PROTEIN/                     # Protein alignments
│
├── EVM/                         # EvidenceModeler output
│
├── FILTER/
│   ├── portcullis/              # Junction filtering
│   ├── mikado/                  # Mikado results
│   └── filter.gff3              # Filter phase output
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
  --out MySpecies_v1.0 --splice-name t --justify 5 --sort
```

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
| Out of memory | Increase `memory` in cluster config for the rule |
| File not found (Singularity) | Ensure paths are within working directory or use `SINGULARITY_BIND` |
| SLURM account error | Use `account` (billing account), not username |
| LFS files not downloaded | Run `git lfs pull` |

### Memory Guidelines

- Recommend **4GB per thread**
- Example: 48 threads = 192g memory
- `ncpus` and `threads` should match

## No HPC?

Deploy a SLURM cluster on Google Cloud: [Cloud Cluster Toolkit](https://docs.cloud.google.com/cluster-toolkit/docs/quickstarts/slurm-cluster)

## Citation

> Sylvan: A comprehensive genome annotation pipeline. *Under review.*

## License

MIT License - see [LICENSE](LICENSE)

## Contact

Issues: https://github.com/plantgenomicslab/Sylvan/issues
