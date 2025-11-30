# Sylvan Genome Annotation

Sylvan is a comprehensive genome annotation pipeline that combines EVM/PASA, GETA, Mikado, and Helixer with semi-supervised filtering for high-quality gene models.

![Sylvan Workflow](https://github.com/plantgenomicslab/Sylvan/blob/main/docs/images/sylvan_workflow.jpg?raw=true)

## Pipeline Workflow

The Sylvan pipeline consists of interconnected modules that process evidence from multiple sources and combine them into a unified gene model.

![Pipeline DAG](https://github.com/plantgenomicslab/Sylvan/blob/main/docs/images/rulegraph.png?raw=true)

The workflow includes:
- **Repeat Masking**: RepeatMasker with species-specific and custom libraries
- **RNA-seq Processing**: STAR alignment, StringTie/PsiCLASS assembly, PASA refinement
- **Protein Homology**: Miniprot, GeneWise, and GMAP alignments
- **Ab Initio Prediction**: Helixer (deep learning) and Augustus (HMM-based)
- **Liftover**: LiftOff for transferring annotations from neighbor species
- **Evidence Combination**: EvidenceModeler (EVM) integration with PASA updates
- **Consistent outputs**: All intermediate/final outputs live under `<repo>/results/` (e.g., `results/EVM`, `results/GETA/RepeatMasker/...`); no working-directory-sensitive symlinks.

## Features

- **Multi-evidence integration**: RNA-seq, protein homology, neighbor species annotations
- **Multiple ab initio predictors**: Helixer, Augustus, GETA
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

## Input Requirements

### Annotation Phase

| Input | Description | Config Field |
|-------|-------------|--------------|
| Genome assembly | FASTA file (`.fa`, `.fasta`, `.fa.gz`, `.fasta.gz`) | `genome` |
| RNA-seq data | Gzipped FASTQ files in a folder | `rna_seq` |
| Protein sequences | FASTA from UniProt, OrthoDB, etc. | `proteins` |
| Neighbor species | GFF3 + genome FASTA files | `liftoff.neighbor_gff`, `liftoff.neighbor_fasta` |
| Repeat library | EDTA output (`.TElib.fa`) | `geta.RM_lib` |
| Singularity image | Path to `sylvan.sif` | `singularity` |
| RexDB (filter) | RepeatExplorer protein DB (retroelement proteins); e.g. Viridiplantae_v4.0.fasta from https://github.com/repeatexplorer/rexdb | `RexDB` |

### Filter Phase (additional)

| Input | Description | Config Field |
|-------|-------------|--------------|
| RexDB | Plant repeat database (e.g. Viridiplantae_v4.0.fasta) | `RexDB` |
| BUSCO lineage | e.g., `eudicots_odb10` | `busco_lin` |

## Configuration

Sylvan uses two separate configuration files:

| File | Purpose |
|------|---------|
| `config_annotate.yml` | **Pipeline options**: input paths, species parameters, tool settings |
| `cluster_annotate.yml` | **SLURM resources**: CPU, memory, partition for each rule |

**`config_annotate.yml`** contains:
- Input file paths (genome, RNA-seq, proteins, neighbor species)
- Species-specific settings (Helixer model, Augustus species)
- Tool parameters (max intron length, EVM weights)
- Output prefix and directories

**`cluster_annotate.yml`** contains:
- SLURM account and partition (in `__default__` section)
- Per-rule resource allocation (CPUs, memory, time limit)
- Log file locations

Example structure:
```yaml
__default__:
  account: your-account      # SLURM billing account
  partition: your-partition  # SLURM partition (queue)
  memory: 4g                 # Default memory
  ncpus: 1                   # Default CPUs
  time: 14-00:00:00          # Max wall time
  output: results/logs/{rule}_{wildcards}.out
  error: results/logs/{rule}_{wildcards}.err

# Override per rule:
pasa:
  ncpus: 64
  threads: 64
  memory: 216g

Output locations:
- `results/` is created under the repo root; all modules write there regardless of current working directory.
- RepeatMasker/RepeatModeler run inside `results/GETA/RepeatMasker/...`, so temporary `.RepeatMaskerCache`/`RM_*` folders also stay under `results/`.
- EVM commands/outputs live in `results/EVM/` (absolute paths in command lists; no `EVM -> results/EVM` symlink needed).
```

This separation allows you to reuse the same pipeline config across different clusters by only changing the cluster config.

Copy and edit the configuration files:

```bash
cp config/config_annotate.yml my_config.yml
cp config/cluster_annotate.yml my_cluster.yml
```

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

## Running the Pipeline

### Annotation Phase

```bash
# Dry run
snakemake -n --snakefile bin/Snakefile_annotate

# Submit to SLURM
sbatch -A [account] -p [partition] -c 1 --mem=1g \
  -J annotate -o annotate.out -e annotate.err \
  --wrap="./bin/annotate.sh"
```

### Filter Phase

```bash
sbatch -A [account] -p [partition] -c 1 --mem=4g \
  -J filter -o filter.out -e filter.err \
  --wrap="./bin/filter.sh"
```

### Useful Commands

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

All outputs are in `results/`:

```
results/
├── complete_draft.gff3      # Final combined annotation
├── AB_INITIO/Helixer/       # Helixer predictions
├── GETA/Augustus/           # Augustus predictions
├── LIFTOVER/LiftOff/        # Neighbor species liftover
├── TRANSCRIPT/PASA/         # PASA assemblies
├── PROTEIN/                 # Protein alignments
├── EVM/                     # EvidenceModeler output
├── FILTER/filter.gff3       # Filtered final output
└── logs/                    # SLURM job logs
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
