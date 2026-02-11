# Sylvan Tutorial: Running with Toy Data

This tutorial walks you through running Sylvan on the included toy dataset (*A. thaliana* chromosome 4). The commands and parameter choices shown below are **specific to this experiment**; when annotating your own genome, adjust the configuration files to match your organism, data, and cluster environment. For a description of the pipeline's architecture and available modules, see the [README](README.md#pipeline-design).

## Table of Contents

- [Prerequisites](#prerequisites)
- [Step 1: Environment Setup](#step-1-environment-setup)
- [Step 2: Download Containers](#step-2-download-containers)
- [Step 3: Prepare Repeat Library (EDTA)](#step-3-prepare-repeat-library-edta)
- [Step 4: Run Annotation Pipeline](#step-4-run-annotation-pipeline)
- [Step 5: Run Filter Pipeline](#step-5-run-filter-pipeline)
- [Step 5b: Alternative Score-based Filter](#step-5b-alternative-score-based-filter)
- [Step 6: Format Output (TidyGFF)](#step-6-format-output-tidygff)
- [Step 7: Cleanup Intermediate Files](#step-7-cleanup-intermediate-files)
- [Advanced Configuration](#advanced-configuration)
- [Monitoring and Debugging](#monitoring-and-debugging)
- [Toy Data Details](#toy-data-details)

---

## Prerequisites

- Linux system with SLURM
- Singularity 3.x+
- Conda/Mamba
- Git LFS

### Host Dependencies

Most bioinformatics tools run inside the Singularity container. The host environment needs:

| Package | Purpose |
|---------|---------|
| Python 3.10+ | Pipeline orchestration and filter scripts |
| Snakemake 7 | Workflow engine |
| pandas | Data manipulation in Snakemake and filter scripts |
| scikit-learn | Random forest classifier (`Filter.py`, `score_filter.py`) |
| NumPy | Numerical operations |
| PyYAML | Config file parsing |
| rich | Logging format (optional) |

Perl (`fillingEndsOfGeneModels.pl`) and R (`filter_distributions.R`) scripts run inside the container and do not need host installation.

## Step 1: Environment Setup

```bash
# Create conda environment
conda create -n sylvan -c conda-forge -c bioconda python=3.11 snakemake=7 -y
conda activate sylvan

# Install Git LFS
git lfs install
```

## Step 2: Download Containers

### Sylvan Container

```bash
singularity pull --arch amd64 sylvan.sif library://wyim/sylvan/sylvan:latest
```

### EDTA Container (for repeat library)

```bash
export SINGULARITY_CACHEDIR=$PWD
singularity pull EDTA.sif docker://quay.io/biocontainers/edta:2.2.0--hdfd78af_1
```

### Clone Repository

```bash
git clone https://github.com/plantgenomicslab/Sylvan.git
cd Sylvan

# Verify LFS files are downloaded (not pointers)
git lfs pull
ls -la toydata/  # Files should be > 200 bytes
```

## Step 3: Prepare Repeat Library (EDTA)

> **Note:** The toy data includes a pre-computed repeat library. This step is for reference or if you need to regenerate it.

```bash
sbatch -c 16 --mem=68g --wrap="singularity exec --cleanenv --env PYTHONNOUSERSITE=1 \
  EDTA.sif EDTA.pl \
  --genome toydata/genome_input/genome.fasta \
  --cds toydata/cds_aa/neighbor.cds \
  --anno 1 --threads 16 --force 1"
```

**Expected runtime:** ~1.5 hours for 18.6 Mb genome

**Output:** `genome.fasta.mod.EDTA.TElib.fa`

### EDTA Benchmark (Toy Data)

| Stage | Duration |
|-------|----------|
| LTR detection | ~3 min |
| SINE detection | ~6 min |
| LINE detection | ~45 min |
| TIR detection | ~8 min |
| Helitron detection | ~10 min |
| Filtering & annotation | ~8 min |
| **Total** | **~1.5 hours** |

## Experiment Setup (Toy Data)

The toy data experiment uses the following configuration choices. These are pre-set in `toydata/config/config_annotate.yml` and `toydata/config/config_filter.yml`:

| Parameter | Choice for this experiment | Alternatives |
|-----------|---------------------------|--------------|
| **Genome** | *A. thaliana* Chr4, 18.6 Mb, 3 segments | Any FASTA assembly |
| **Repeat masking** | RepeatMasker with `Embryophyta` library + custom EDTA library | Any RepeatMasker species; RepeatModeler for de novo |
| **RNA-seq** | 12 paired-end samples (leaf, rosette, whole plant) | Any number of paired-end samples |
| **RNA-seq alignment** | STAR pathway with StringTie + PsiCLASS | HiSat2 is also available |
| **Protein homology** | Combined neighbor-species proteins (`neighbor.aa`) | UniProt, OrthoDB, or any protein FASTA |
| **Neighbor species** | 3 species: *A. lyrata*, *C. rubella*, *C. hirsuta* | One or more annotated relatives |
| **Ab initio (Helixer)** | `land_plant` model, subsequence length 64152 | `vertebrate`, `invertebrate`, `fungi` |
| **Ab initio (Augustus)** | Start from `arabidopsis` model, train on target data | Any existing Augustus species; or skip training |
| **EVM weights** | PASA=10, Augustus=7, GETA=5, Helixer=3 | Adjust to evidence quality |
| **Filter cutoffs** | TPM=3, coverage=0.5, BLAST identity=0.6 | Adjust per organism |
| **BUSCO lineage** | `eudicots_odb10` | Any BUSCO lineage |

These choices are experiment-specific. The pipeline supports all the alternatives listed above — see the [README](README.md#pipeline-design) for the full set of configurable modules.

---

## Step 4: Run Annotation Pipeline

> **Note:** The commands below use the pre-configured toy data settings. For your own data, replace `toydata/config/config_annotate.yml` with your config file path.

### Environment Variables

The pipeline uses several environment variables for configuration:

| Variable | Description | Example |
|----------|-------------|---------|
| `SYLVAN_CONFIG` | Path to pipeline config file (also used as `--cluster-config` for SLURM settings) | `toydata/config/config_annotate.yml` |
| `SYLVAN_FILTER_CONFIG` | Path to filter config file | `toydata/config/config_filter.yml` |
| `SYLVAN_RESULTS_DIR` | Override results output directory (default: `$(pwd)/results/`). Useful on HPC systems with separate storage. | `/scratch/$USER/sylvan_results` |
| `TMPDIR` | Temporary directory for intermediate files | `$(pwd)/results/TMP` |
| `SLURM_TMPDIR` | SLURM job temporary directory (should match `TMPDIR`) | `$TMPDIR` |

**Why set `TMPDIR`?**

Setting `TMPDIR` explicitly is **critical** on many HPC systems:

1. **Memory-backed `/tmp` (tmpfs)**: Some HPC nodes have no local disk storage. The default `/tmp` is mounted as `tmpfs`, which stores files **in RAM**. Large temporary files from tools like STAR, RepeatMasker, or Augustus can quickly exhaust memory and crash your jobs with cryptic "out of memory", "no space left on device", or even **segmentation fault** errors—since the OS may kill processes or corrupt memory when tmpfs fills up.

2. **Quota limits**: Shared `/tmp` partitions often have strict per-user quotas (e.g., 1-10 GB). Genome annotation tools easily exceed this.

3. **Job isolation**: When `TMPDIR` points to your project directory, temp files persist after job completion for debugging. Cleanup is also straightforward with `rm -rf results/TMP/*`.

4. **Singularity compatibility**: Containers inherit `TMPDIR` from the host. Setting it to a bound path ensures temp files are written to accessible storage.

> **Tip**: Check if your cluster uses tmpfs with `df -h /tmp`. If it shows `tmpfs` as the filesystem type, you **must** set `TMPDIR` to avoid memory issues.

### Complete Environment Setup

Copy this block before running the pipeline:

```bash
# Required: Pipeline configuration
export SYLVAN_CONFIG="toydata/config/config_annotate.yml"

# Required: Temp directory (create if not exists)
mkdir -p results/TMP
export TMPDIR="$(pwd)/results/TMP"
export SLURM_TMPDIR="$TMPDIR"

# Optional: Bind additional paths for Singularity
# export SINGULARITY_BIND="/scratch,/data"

# Optional: Increase open file limit (some tools need this)
ulimit -n 65535 2>/dev/null || true
```

### Additional Singularity Variables

| Variable | Description | When to Use |
|----------|-------------|-------------|
| `SINGULARITY_BIND` | Bind additional host paths into container | When input files are outside working directory |
| `SINGULARITY_CACHEDIR` | Location for Singularity cache | When home directory has quota limits |
| `SINGULARITY_TMPDIR` | Singularity's internal temp directory | Should match `TMPDIR` |

**When do you need `SINGULARITY_BIND`?**

Singularity automatically binds your current working directory, `$HOME`, and `/tmp`. However, you need explicit binding when:

- Input files (genome, RNA-seq, proteins) are on a **separate filesystem** (e.g., `/scratch`, `/project`, `/data`)
- Your **home directory has quota limits** and you store data elsewhere
- Using **shared lab storage** mounted at non-standard paths
- Config file references **absolute paths** outside the working directory

```bash
# Common scenarios:

# 1. Data on scratch space
export SINGULARITY_BIND="/scratch/$USER"

# 2. Multiple paths (comma-separated)
export SINGULARITY_BIND="/scratch,/project/shared_data,/data/genomes"

# 3. Bind with different container path (host:container)
export SINGULARITY_BIND="/long/path/on/host:/data"

# 4. Read-only binding (for shared reference data)
export SINGULARITY_BIND="/shared/databases:/databases:ro"
```

**Diagnosing bind issues:**
```bash
# Error: "file not found" inside container but exists on host
# → The path isn't bound. Add it to SINGULARITY_BIND

# Test if path is accessible inside container:
singularity exec sylvan.sif ls /your/data/path

# See what's currently bound:
singularity exec sylvan.sif cat /proc/mounts | grep -E "scratch|project|data"
```

### Verifying Your Environment

Before submitting jobs, verify your setup:

```bash
# Check TMPDIR is on real disk (not tmpfs)
df -h $TMPDIR

# Verify Singularity can access paths
singularity exec sylvan.sif ls $TMPDIR

# Test config file is readable
cat $SYLVAN_CONFIG | head -5
```

### Debugging Environment Issues

When jobs fail unexpectedly, environment variables are often the culprit. Use these techniques to diagnose:

**1. Print all relevant variables:**
```bash
# Add to your script or run interactively
echo "=== Environment Check ==="
echo "SYLVAN_CONFIG: $SYLVAN_CONFIG"
echo "TMPDIR: $TMPDIR"
echo "SLURM_TMPDIR: $SLURM_TMPDIR"
echo "SINGULARITY_BIND: $SINGULARITY_BIND"
echo "PWD: $PWD"
df -h $TMPDIR
```

**2. Check what SLURM jobs actually see:**
```bash
# Submit a diagnostic job
sbatch --wrap='env | grep -E "TMPDIR|SINGULARITY|SYLVAN" && df -h /tmp $TMPDIR'
```

**3. Common environment-related errors:**

| Error Message | Likely Cause | Solution |
|---------------|--------------|----------|
| `No space left on device` | TMPDIR on tmpfs or quota exceeded | Set `TMPDIR` to project storage |
| `Segmentation fault` | Memory exhausted (tmpfs full) | Set `TMPDIR` to disk-backed storage |
| `file not found` (in container) | Path not bound in Singularity | Add path to `SINGULARITY_BIND` |
| `Permission denied` | Singularity can't write to TMPDIR | Check directory permissions, ensure path is bound |
| `cannot create temp file` | TMPDIR doesn't exist or not writable | Run `mkdir -p $TMPDIR && touch $TMPDIR/test` |

**4. Interactive debugging inside container:**
```bash
# Start interactive shell in container with same bindings
singularity shell --cleanenv sylvan.sif

# Inside container, verify paths exist
ls -la $TMPDIR
ls -la /path/to/your/data
```

**5. Check if variables survive into SLURM jobs:**

SLURM doesn't always pass environment variables. Ensure your submit script exports them:
```bash
# In your sbatch script or wrapper
#SBATCH --export=ALL    # Pass all environment variables

# Or explicitly export in the script:
export TMPDIR="$(pwd)/results/TMP"
export SLURM_TMPDIR="$TMPDIR"
```

### Dry Run (Recommended)

Always do a dry run first to verify configuration:

```bash
export SYLVAN_CONFIG="toydata/config/config_annotate.yml"
snakemake -n --snakefile bin/Snakefile_annotate
```

### Submit to SLURM

```bash
sbatch -A [account] -p [partition] -c 1 --mem=1g \
  -J annotate -o annotate.out -e annotate.err \
  --wrap="./bin/annotate_toydata.sh"
```

**Output locations**  
- All intermediate/final results are written under the repo root `results/`.  
- RepeatMasker/RepeatModeler run inside `results/GETA/RepeatMasker/...`, so `.RepeatMaskerCache` and `RM_*` temp folders also stay there.  
- EVM commands and outputs live in `results/EVM/`; no `EVM -> results/EVM` symlink is needed.
- For the filter pipeline, set `RexDB` to a RepeatExplorer protein DB (e.g., Viridiplantae_v4.0.fasta from https://github.com/repeatexplorer/rexdb). You can download directly via:  
  `wget -O toydata/misc/Viridiplantae_v4.0.fasta https://raw.githubusercontent.com/repeatexplorer/rexdb/refs/heads/main/Viridiplantae_v4.0.fasta`

### Expected Runtime (Toy Data)

| Stage | Time |
|-------|------|
| RepeatMasking | 15-30 min |
| RNA-seq alignment | 30-60 min |
| Transcript assembly | 20-40 min |
| Homology search | 30-60 min |
| Augustus training | 1-2 hours |
| Gene model combination | 30-60 min |
| EvidenceModeler | 30-60 min |
| **Total** | **4-8 hours** |

### Force Rerun

```bash
# Rerun all jobs
./bin/annotate_toydata.sh --forceall

# Rerun specific rule
./bin/annotate_toydata.sh --forcerun helixer

# Rerun incomplete jobs
./bin/rerun-incomplete_toydata.sh
```

## Step 5: Run Filter Pipeline

> **Note:** The filter cutoff thresholds and parameters below are configured for the toy data experiment. Adjust `config_filter.yml` for your organism.

After annotation completes:

```bash
# Set filter config
export SYLVAN_FILTER_CONFIG="toydata/config/config_filter.yml"

# Dry run
snakemake -n --snakefile bin/Snakefile_filter

# Submit to SLURM
sbatch -A [account] -p [partition] -c 1 --mem=4g \
  -J filter -o filter.out -e filter.err \
  --wrap="./bin/filter_toydata.sh"
```

### Filter Pipeline Steps

The filter phase runs these steps in order:

1. **Extract sequences**: CDS and peptide extraction from annotated GFF (`gff3_file_to_proteins.pl`)
2. **PfamScan**: Identify conserved protein domains using Pfam-A HMMs
3. **RSEM**: Quantify transcript expression (TPM) from re-aligned RNA-seq reads
4. **BLASTp**: Homolog similarity (protein DB) and repeat similarity (RexDB) — parallelized across 20 split peptide files
5. **Ab initio coverage**: Calculate overlap with Augustus, Helixer, and RepeatMasker predictions
6. **lncDC**: XGBoost-based classification of long non-coding RNAs using plant-specific pre-trained models
7. **BUSCO**: Identify conserved gene models (monitoring only, not used as a filter feature)
8. **Semi-supervised random forest**: Iterative training with configurable thresholds

### Filter Cutoff Thresholds

Thresholds in `config_filter.yml` under `Cutoff` control the initial heuristic gene selection:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `tpm` | TPM threshold | 3 |
| `rsem_cov` | RNA-seq coverage | 0.5 |
| `blast_pident` / `blast_qcovs` | BLASTp identity/coverage | 0.6 / 0.6 |
| `rex_pident` / `rex_qcovs` | RexDB identity/coverage | 0.6 / 0.6 |
| `helixer_cov` / `augustus_cov` | Ab initio overlap | 0.8 / 0.8 |
| `repeat_cov` | Repeat overlap | 0.5 |

### Filter Random Forest Parameters

Set via `Filter.py` arguments (configured in `Snakefile_filter`):

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--trees` | Number of trees in the random forest | 100 |
| `--predictors` | Number of predictors for tree splitting | 5 |
| `--recycle` | Predicted accuracy required to add observation to next iteration | 0.95 |
| `--max-iter` | Maximum re-training iterations | 10 |
| `--seed` | Random seed for reproducibility | 123 |

The filter typically converges within 3-5 iterations. Each iteration adds high-confidence predictions back into the training set and retrains. The process stops when no new observations exceed the `--recycle` threshold or `--max-iter` is reached.

### Important: `chrom_regex`

The `chrom_regex` field in `config_filter.yml` is required for proper chromosome identification. It must match your genome's chromosome naming convention:

```yaml
# Common patterns:
chrom_regex: (^Chr)|(^chr)|(^LG)|(^Ch)|(^\d)
```

### Output Files

| File | Description |
|------|-------------|
| `results/FILTER/filter.gff3` | Kept gene models |
| `results/FILTER/discard.gff3` | Discarded gene models |
| `results/FILTER/data.tsv` | Feature matrix used by random forest |
| `results/FILTER/keep_data.tsv` | Evidence data for kept genes |
| `results/FILTER/discard_data.tsv` | Evidence data for discarded genes |
| `results/FILTER/{prefix}.cdna` | Transcript sequences |
| `results/FILTER/{prefix}.pep` | Peptide sequences |
| `results/FILTER/lncrna_predict.csv` | lncDC predictions |

## Step 5b: Alternative Score-based Filter

An alternative scoring pipeline (`Snakefile_filter_score`) is available. Instead of the iterative semi-supervised approach, it uses logistic regression and random forest scoring with pseudo-labels derived from evidence data:

- **Positive pseudo-labels**: genes with Pfam hits or homolog BLAST hits
- **Negative pseudo-labels**: genes with RexDB (repeat) hits
- **Threshold selection**: maximizes F1 score on pseudo-labels via precision-recall curve

```bash
export SYLVAN_FILTER_CONFIG="toydata/config/config_filter.yml"
./bin/filter_score_toydata.sh
```

**Output:**
- `results/FILTER/scores.csv` — Per-gene scores, features, and predictions
- `results/FILTER/scores.metrics.txt` — AUC, precision-recall, F1, and chosen thresholds

This approach is simpler and faster but may be less accurate than the semi-supervised method for organisms with limited prior annotation data.

## Step 6: Format Output (TidyGFF)

```bash
singularity exec sylvan.sif python bin/TidyGFF.py \
  Ath4 results/FILTER/filter.gff3 \
  --out Ath4_v1.0 \
  --splice-name t \
  --justify 5 \
  --sort \
  --chrom-regex "^Chr" \
  --source Sylvan
```

## Step 7: Cleanup Intermediate Files

After **both** annotation and filter phases have completed successfully, run the cleanup script to remove intermediate files:

```bash
./bin/cleanup.sh
```

**What it removes:**
- Untracked Snakemake workflow files from the annotation phase
- Trimmed RNA-seq files (`GETA/fastp/single/*.fq.gz`, `GETA/fastp/paired/*.fq.gz`)
- RepeatMasker cache directories (`.RepeatMaskerCache/`)
- Empty directories
- Temporary PASA and HMM files

**What it preserves:**
- Final outputs (`complete_draft.gff3`, `FILTER/`)
- Log files (`results/logs/`)
- Configuration files, Singularity images, and pipeline scripts
- LncDC database files

> **Warning:** Only run this after you have verified both pipeline phases completed successfully. The annotation phase intermediate files cannot be regenerated without re-running the full annotation pipeline.

---

## Advanced Configuration

### EVM Weights (`evm_weights.txt`)

Controls how EvidenceModeler prioritizes evidence sources. Higher weights give more influence:

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

- **PASA transcripts** (weight 10) have the highest weight as direct transcript evidence
- **Augustus** (weight 7) is weighted higher than **Helixer** (weight 3) as the primary ab initio predictor
- **GETA** combined models (weight 5) reflect the integrated GETA pipeline evidence
- Adjust weights based on evidence quality for your organism

### EVM Partitioning (`num_evm_files`)

The `num_evm_files` parameter in `config_annotate.yml` controls how many parallel EVM partitions to create. The genome is split into overlapping segments (1 Mb segments with 20 kb overlap):

- **Higher values** = more parallel SLURM jobs = faster wall-clock time, but more job scheduling overhead
- **Lower values** = fewer jobs = less cluster burden, but slower
- Default: `126` (works well for genomes up to ~500 Mb)
- For very large genomes (>1 Gb), consider increasing to 200-500

### Mikado Scoring (`config/plant.yaml`)

The `plant.yaml` file controls Mikado transcript selection parameters. Default values are tuned for plants with intron sizes typical of most angiosperms. Key parameters:

- **requirements**: Minimum CDS fraction, exon count, intron length constraints
- **not_fragmentary**: Criteria for rejecting incomplete models
- **scoring**: Weights for transcript ranking (blast score, CDS length, UTR length, etc.)

Most users should not need to modify this file. For organisms with unusual intron sizes or gene structures, consult the [Mikado documentation](https://mikado.readthedocs.io/).

### Augustus Training Options

| Config Parameter | Description |
|-----------------|-------------|
| `augustus_species` | Species name for Augustus training (e.g., `arabidopsis`) |
| `augustus_start_from` | Start training from an existing species model (faster convergence for closely related species) |
| `use_augustus` | Use a pre-trained Augustus model without re-training (set to species name; `placeholder` = train fresh) |

- If your organism is close to a well-annotated species, set `augustus_start_from` to that species for faster, more accurate training
- If Augustus training fails (requires ~500 training genes minimum), use `use_augustus` with a close species

### Helixer GPU Configuration

Helixer runs ~10x faster with GPU acceleration. Configure a separate GPU partition in the per-rule overrides section of `config_annotate.yml`:

```yaml
helixer:
  ncpus: 4
  memory: 32g
  account: your-gpu-account
  partition: your-gpu-partition
```

### RNA-seq File Naming

The pipeline supports paired-end reads with two naming conventions:
- `{sample}_1.fastq.gz` / `{sample}_2.fastq.gz`
- `{sample}_R1.fastq.gz` / `{sample}_R2.fastq.gz`

All paired-end FASTQ files should be placed in the directory specified by `rna_seq` in the config. Sample names are automatically detected from file names.

### Annotation Pipeline Architecture

The annotation Snakefile supports two parallel RNA-seq alignment pathways:

1. **STAR pathway**: `STAR_paired` → `psiClass_STAR` / `stringtie_STAR`
2. **HiSat2 pathway**: `HiSat2_PAIRED` → `psiClass_HiSat` / `stringtie_HiSat`

Both pathways feed into transcript assembly and the GETA pipeline. The choice is determined by the config and available resources. STAR is generally preferred for accuracy; HiSat2 uses less memory.

The protein homology pipeline is sequential: **Miniprot** (fast protein-to-genome alignment) → **miniprot2genewise** (converts to gene regions) → **GeneWise** (refined gene structure prediction on identified regions).

### Snakemake Job Grouping

The `annotate.sh` script uses Snakemake job grouping for the `Sam2Transfrag` rule:

```bash
--groups Sam2Transfrag=group0 --group-components group0=100
```

This groups up to 100 `Sam2Transfrag` jobs into a single SLURM submission, reducing cluster scheduling overhead for this highly parallelized step. Adjust the group-components value if you experience SLURM scheduling issues or want different parallelism.

---

## Monitoring and Debugging

### Check Job Status

```bash
squeue -u $USER
```

### View Logs

```bash
# Snakemake log
tail -f .snakemake/log/*.snakemake.log

# Find recent error logs
ls -lt results/logs/*.err | head -10

# Search for errors
grep -l 'Error\|Traceback' results/logs/*.err

# View specific log (pattern: {rule}_{wildcards}.err)
cat results/logs/liftoff_.err
cat results/logs/geneRegion2Genewise_seqid=group17400.err
```

### Common Issues

| Issue | Solution |
|-------|----------|
| Out of memory | Increase `memory` in the per-rule overrides section of `config_annotate.yml` for the failing rule |
| `No space left on device` | `TMPDIR` is on tmpfs or quota exceeded — set `TMPDIR` to project storage |
| `Segmentation fault` | Often caused by tmpfs exhaustion — set `TMPDIR` to disk-backed storage |
| LFS files are pointers | Run `git lfs pull`; verify with `ls -la toydata/` (files should be > 200 bytes) |
| Singularity bind error | Ensure paths are within working directory or use `SINGULARITY_BIND` |
| Augustus training fails | Needs ~500 training genes minimum; use `augustus_start_from` with a close species or `use_augustus` to skip training |
| Job timeout | Increase `time` in the per-rule overrides section of `config_annotate.yml` for the rule |
| Variables not in SLURM job | Add `#SBATCH --export=ALL` or explicitly export in submit script |
| Filter `chrom_regex` error | Ensure `chrom_regex` in `config_filter.yml` matches your chromosome naming convention |

---

## Toy Data Details

### Overview

The toy dataset contains *Arabidopsis thaliana* chromosome 4 split into 3 segments (~18.6 Mb total).

### Directory Structure

```
toydata/
├── config/                      # Configuration files
│   ├── config_annotate.yml      # Annotation pipeline + SLURM resource config (pre-configured)
│   ├── config_filter.yml        # Filter pipeline config
│   ├── evm_weights.txt          # EVM evidence weights
│   └── plant.yaml               # Mikado transcript scoring parameters
├── genome_input/
│   └── genome.fasta.gz          # A. thaliana Chr4 (3 parts)
├── RNASeq/                      # 12 paired-end RNA-seq samples
│   ├── sub_SRR1019221_1.fastq.gz
│   └── ...
├── neighbor_genome/             # Neighbor species genomes
│   ├── aly4.fasta               # A. lyrata
│   ├── ath4.fasta               # A. thaliana
│   ├── chi4.fasta               # C. hirsuta
│   └── cru4.fasta               # C. rubella
├── neighbor_gff3/               # Neighbor annotations
│   ├── aly4.gff3
│   └── ...
├── cds_aa/
│   ├── neighbor.cds             # Combined CDS for EDTA
│   └── neighbor.aa              # Proteins for homology
├── misc/
│   └── Viridiplantae_v4.0.fasta # RexDB plant repeat protein database
└── EDTA/                        # Pre-computed repeat library
    └── genome.fasta.mod.EDTA.TElib.fa
```

### Genome Statistics

| Segment | Length | GC (%) |
|---------|--------|--------|
| Chr4_1 | 6,195,060 bp | 36.66 |
| Chr4_2 | 6,195,060 bp | 35.24 |
| Chr4_3 | 6,195,018 bp | 36.69 |
| **Total** | **18,585,138 bp** | 36.20 |

### TAIR10 Reference Annotation (Chr4)

For comparison, the official TAIR10 annotation of *Arabidopsis thaliana* (Col-0) chromosome 4 contains:

| Feature Type | Count |
|--------------|-------|
| Protein-coding genes | 4,124 |
| pre-tRNA genes | 79 |
| rRNA genes | 0 |
| snRNA genes | 0 |
| snoRNA genes | 11 |
| miRNA genes | 28 |
| Other RNA genes | 62 |
| Pseudogenes | 121 |
| Transposable element (TE) genes | 711 |
| **Total annotated loci** | **5,410** |

> **Note:** TAIR10 is the current reference annotation standard for *A. thaliana*. The protein-coding gene count increased from 3,744 (original Chr4 paper) to 4,124 as annotation methods improved.
>
> Source: [Phoenix Bioinformatics - Genome Annotation at TAIR](https://phoenixbioinformatics.atlassian.net/wiki/spaces/COM/pages/42216279/Genome+Annotation+at+TAIR)

### Sylvan Pipeline Results (Toy Data)

Running the Sylvan pipeline on the Chr4 toy dataset produces the following results:

**Annotation Phase (complete_draft.gff3):**

| Metric | Count |
|--------|-------|
| Total genes | 5,720 |
| Total mRNA | 5,800 |

**Filter Phase (filter.gff3):**

| Metric | Count |
|--------|-------|
| Genes kept | 3,756 |
| mRNA kept | 3,834 |
| Genes discarded | 1,964 |
| mRNA discarded | 1,966 |

**Output files:**
- `results/FILTER/filter.gff3` - Kept gene models
- `results/FILTER/discard.gff3` - Discarded gene models
- `results/FILTER/data.tsv` - Feature matrix used by random forest (input to feature importance analysis)
- `results/FILTER/keep_data.tsv` - Evidence data for kept genes
- `results/FILTER/discard_data.tsv` - Evidence data for discarded genes
- `results/FILTER/{prefix}.cdna` - Extracted transcript sequences
- `results/FILTER/{prefix}.pep` - Extracted peptide sequences
- `results/FILTER/lncrna_predict.csv` - lncDC long non-coding RNA predictions
- `results/complete_draft.gff3.map` - ID mapping between original and new IDs

> **Comparison:** The 3,756 kept genes represents ~91% of TAIR10's 4,124 protein-coding genes on Chr4. The higher initial count (5,720) includes transposable elements (TAIR10 has 711 TE genes) and low-confidence predictions that are filtered out.

### Neighbor Species

| Code | Species | Common Name |
|------|---------|-------------|
| aly4 | *Arabidopsis lyrata* | Lyrate rockcress |
| cru4 | *Capsella rubella* | Pink shepherd's purse |
| chi4 | *Cardamine hirsuta* | Hairy bittercress |

### RNA-seq Samples

| SRA | Tissue | Size |
|-----|--------|------|
| SRR1019221 | Leaf (14-day) | 4.6 Gb |
| SRR1105822 | Rosette (19-day) | 3.1 Gb |
| SRR1105823 | Rosette (19-day) | 4.5 Gb |
| SRR1106559 | Rosette (19-day) | 3.6 Gb |
| SRR446027 | Whole plant | 2.5 Gb |
| SRR446028 | Whole plant | 5.4 Gb |
| SRR446033 | Whole plant | 5.3 Gb |
| SRR446034 | Whole plant | 5.4 Gb |
| SRR446039 | Whole plant | 2.5 Gb |
| SRR446040 | Whole plant | 6.3 Gb |
| SRR764885 | Leaf (4-week) | 4.6 Gb |
| SRR934391 | Whole plant | 4.0 Gb |

### How the Toy Data Was Created

**1. Genome segmentation**
```bash
seqkit split2 -p 3 Chr4.fasta
```

**2. Neighbor CDS extraction**

Syntenic regions were identified using [MCscan/jcvi](https://github.com/tanghaibao/jcvi):
```bash
python -m jcvi.compara.catalog ortholog ath aly --no_strip_names
python -m jcvi.compara.synteny mcscan ath.bed ath.aly.lifted.anchors --iter=1
```

**3. RNA-seq subsetting**

Reads mapping to Chr4 were extracted:
```bash
STAR --genomeDir star_index --readFilesIn sample_1.fq.gz sample_2.fq.gz
samtools view -b -F 4 Aligned.bam | samtools fastq -1 out_1.fq -2 out_2.fq -
```

### Test Environment

The toy data was tested on:

| Specification | Value |
|---------------|-------|
| Nodes | 4 |
| Total CPUs | 256 |
| CPU | Intel Xeon E5-2683 v4 @ 2.10GHz |
| Cores per node | 64 (2 sockets × 16 cores × 2 threads) |
| Memory per node | 256 GB |
| Storage | GPFS |

### Runtime Statistics

The following summarizes the runtime distribution across all pipeline rules when running on the toy dataset.

**Key observations:**
- **Most time-consuming steps**: `aggregate_CombineGeneModels` (~50,000s), `geneRegion2Genewise` (~1,000s), and `Sam2Transfrag` (~100-200s) are the bottlenecks
- **Fast steps**: Most preprocessing and formatting rules complete in under 10 seconds
- **Parallelizable rules**: Rules like `geneRegion2Genewise`, `gmapExon`, and `STAR_paired` run as multiple parallel jobs (shown as multiple dots), significantly reducing wall-clock time
- **GPU-accelerated**: `helixer` benefits from GPU acceleration when available

**Runtime variability:**

Actual runtime will vary significantly depending on your hardware and cluster configuration:

| Factor | Impact |
|--------|--------|
| **CPU speed** | Faster clock speeds reduce single-threaded bottlenecks |
| **Available nodes** | More nodes = more parallel jobs = faster wall-clock time |
| **Memory per node** | Insufficient memory causes job failures or swapping |
| **Storage I/O** | GPFS/Lustre faster than NFS; SSD faster than HDD |
| **Queue wait time** | Busy clusters add significant delays between jobs |
| **GPU availability** | Helixer runs ~10x faster with GPU acceleration |

With the test environment above (4 nodes, 256 CPUs, 256 GB/node), the toy dataset completes in **4-8 hours** wall-clock time. On smaller clusters or shared resources, expect longer runtimes.

---

## Utility Scripts

### Generate Cluster Config

`bin/generate_cluster_from_config.py` extracts the SLURM-relevant sections (`__default__` and per-rule overrides) from `config_annotate.yml` into a standalone cluster-only YAML. This is optional — `config_annotate.yml` already serves as both pipeline config and `--cluster-config`, so you only need this script if you want a minimal, separate cluster file.

```bash
# Extract standalone cluster YAML from production config
python bin/generate_cluster_from_config.py \
  --config config/config_annotate.yml \
  --out config/cluster_annotate.yml \
  --account your-account --partition your-partition

# Extract from toydata config
python bin/generate_cluster_from_config.py \
  --config toydata/config/config_annotate.yml \
  --out toydata/config/cluster_annotate.yml \
  --account your-account --partition your-partition
```

### Key Python Scripts

| Script | Purpose |
|--------|---------|
| `Filter.py` | Semi-supervised random forest gene model classification |
| `score_filter.py` | Alternative logistic regression + RF scoring pipeline |
| `filter_feature_importance.py` | Leave-one-feature-out ablation study for filter |
| `TidyGFF.py` | Reformat GFF for public distribution (renumber IDs, sort, validate) |
| `CombineDuck.py` | DuckDB-based gene annotation database for model selection |
| `combine_genemodel.py` | Merge Augustus, TransFrag, and GeneWise gene models |
| `Pick_Primaries.py` | Select primary transcript per gene locus |
| `repeat_gene_removal.py` | Remove gene models overlapping with transposable elements |
| `gff_to_evm.py` | Convert GFF3 to EVM input format |
| `splitEVMCommands.py` | Partition EVM commands for parallel execution |
| `clusterGeneWiseRegions.py` | Cluster GeneWise alignment regions |
| `miniprot2Genewise.py` | Convert Miniprot output to GeneWise format |
| `splitBam.py` | Split BAM files for parallel processing |
| `MonitorFilter.py` | Visualize filter iteration progress (matplotlib) |

---

## Getting Help

- Issues: https://github.com/plantgenomicslab/Sylvan/issues
- See also: [README.md](README.md) for configuration reference

## Feature Importance Analysis

After finishing the filter phase you will have `FILTER/data.tsv` (the feature
matrix used by `Filter.py`) and a BUSCO run directory such as
`results/busco/eudicots_odb10`. Reviewers often ask for a feature ablation
study, so we provide an automated helper:

```bash
python bin/filter_feature_importance.py FILTER/data.tsv results/busco/<lineage>/full_table.tsv \
  --output-table FILTER/feature_importance.tsv
```

- **What is the BUSCO full table?** Every BUSCO run writes a
  `full_table.tsv` inside its lineage-specific run folder. Each non-Missing
  BUSCO row lists the BUSCO ID, status (Complete/Duplicated/Fragmented), and the
  transcript/gene ID it matched. The feature-importance script reuses this file
  to count how many BUSCOs remain in the “keep” set during each iteration—no new
  BUSCO analysis is required.
- **Outputs**: `FILTER/feature_importance.tsv` (table) plus
  `FILTER/feature_importance.json` (machine-readable). Both include the baseline
  run (all features) and each leave-one-feature-out run, along with final
  out-of-bag (OOB) error, BUSCO counts, and iteration counts.
- **Optional flags**:
  - `--features TPM COVERAGE PFAM ...` restricts the analysis to specific
    columns from `FILTER/data.tsv`.
  - `--ignore TPM_missing singleExon` removes metadata columns so the script
    automatically uses every other feature column.

Workflow summary:

1. Run `Filter.py` as usual to create `FILTER/data.tsv`.
2. Identify the BUSCO `full_table.tsv` path you already used for filter
   monitoring (e.g., `results/busco/eudicots_odb10/full_table.tsv`).
3. Execute the command above. Inspect `FILTER/feature_importance.tsv` to see how
   dropping each feature affects OOB error (positive delta ⇒ feature is
   important).
4. Incorporate the results (table/plot) into your manuscript or reviewer
   response.
