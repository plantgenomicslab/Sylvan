# Sylvan Tutorial: Running with Toy Data

This tutorial walks you through running Sylvan on the included toy dataset (*A. thaliana* chromosome 4). The commands and parameter choices shown below are **specific to this experiment**; when annotating your own genome, adjust the configuration files to match your organism, data, and cluster environment. For a description of the pipeline's architecture and available modules, see the [README](README.md#pipeline-design).

## Table of Contents

- [Prerequisites](#prerequisites)
- [Toy Data Experiment Overview](#toy-data-experiment-overview)
- [Step 1: Environment Setup](#step-1-environment-setup)
- [Step 2: Download Containers](#step-2-download-containers)
- [Step 3: Prepare Repeat Library (EDTA)](#step-3-prepare-repeat-library-edta)
- [Step 4: Run Annotation Pipeline](#step-4-run-annotation-pipeline)
- [Step 5: Run Filter Pipeline](#step-5-run-filter-pipeline)
- [Step 5b: Alternative Score-based Filter](#step-5b-alternative-score-based-filter)
- [Step 5c: Feature Importance Analysis](#step-5c-feature-importance-analysis)
- [Step 5d: Benchmark (BUSCO + OMArk)](#step-5d-benchmark-busco--omark)
- [Step 6: Format Output (TidyGFF)](#step-6-format-output-tidygff)
- [Step 7: Cleanup Intermediate Files](#step-7-cleanup-intermediate-files)
- [Advanced Configuration](#advanced-configuration)
- [Monitoring and Debugging](#monitoring-and-debugging)
- [Toy Data Details](#toy-data-details)
- [Utility Scripts](#utility-scripts)
- [Runtime Benchmark](#runtime-benchmark)
- [Local Execution (without SLURM)](#local-execution-without-slurm)
- [Getting Help](#getting-help)

---

## Prerequisites

- Linux system (SLURM optional — see [Local Execution](#local-execution-without-slurm))
- Singularity 3.x+
- Conda/Mamba
- Git LFS

For host dependency details, see [README — Dependencies](README.md#dependencies).

---

## Toy Data Experiment Overview

The toy data experiment uses the following configuration choices. These are pre-set in `toydata/config/config_annotate.yml` and `toydata/config/config_filter.yml`:

| Parameter | Choice for this experiment | Alternatives |
|-----------|---------------------------|--------------|
| **Genome** | *A. thaliana* Chr4, 18.6 Mb, 3 segments | Any FASTA assembly (`.fa`, `.fasta`, `.fna`, `.fas`, `.fsa`, `.seq`) |
| **Repeat masking** | RepeatMasker with `Embryophyta` library + custom EDTA library | Any RepeatMasker species; RepeatModeler for de novo |
| **RNA-seq** | 12 paired-end samples (leaf, rosette, whole plant) | Any number of paired-end samples |
| **RNA-seq alignment** | STAR pathway with StringTie + PsiCLASS | HiSat2 is also available |
| **Protein homology** | Combined neighbor-species proteins (`neighbor.aa`) | UniProt, OrthoDB, or any protein FASTA |
| **Neighbor species** | 3 species: *A. lyrata*, *C. rubella*, *C. hirsuta* | One or more annotated relatives |
| **Ab initio (Helixer)** | `land_plant` model, subsequence length 64152 | `vertebrate`, `fungi` |
| **Ab initio (Augustus)** | Start from `arabidopsis` model, train on target data | Any existing Augustus species; or skip training |
| **EVM weights** | PASA=10, Augustus=7, GETA=5, Helixer=3 | Adjust to evidence quality |
| **Filter cutoffs** | TPM=3, coverage=0.5, BLAST identity=0.6 | Adjust per organism |
| **BUSCO lineage** | `eudicots_odb10` | Any BUSCO lineage |

These choices are experiment-specific. The pipeline supports all the alternatives listed above — see the [README](README.md#pipeline-design) for the full set of configurable modules.

---

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
# latest = v4 (GPU-capable TensorFlow, ~16 GB)
singularity pull --arch amd64 sylvan.sif library://wyim/sylvan/sylvan:latest
# Or v3 (CPU-only TensorFlow, ~8 GB):
# singularity pull --arch amd64 sylvan.sif library://wyim/sylvan/sylvan:v3
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

---

## Step 4: Run Annotation Pipeline

> **Note:** The commands below use the pre-configured toy data settings. For your own data, replace `toydata/config/config_annotate.yml` with your config file path.

### Environment Variables

The pipeline uses several environment variables for configuration:

| Variable | Description | Example |
|----------|-------------|---------|
| `SYLVAN_CONFIG` | Path to annotate pipeline config (inputs, tool params, threads) | `toydata/config/config_annotate.yml` |
| `SYLVAN_CLUSTER_CONFIG` | Path to annotate cluster config (SLURM account, partition, resources). Defaults to `SYLVAN_CONFIG` if not set (single-file mode) | `toydata/config/cluster_annotate.yml` |
| `SYLVAN_FILTER_CONFIG` | Path to filter pipeline config | `toydata/config/config_filter.yml` |
| `SYLVAN_FILTER_CLUSTER_CONFIG` | Path to filter cluster config. Defaults to `SYLVAN_FILTER_CONFIG` if not set | `toydata/config/cluster_filter.yml` |
| `SYLVAN_RESULTS_DIR` | Override results output directory (default: `$(pwd)/results/`). Useful on HPC systems with separate storage. | `/scratch/$USER/sylvan_results` |
| `TMPDIR` | Temporary directory for intermediate files | `$(pwd)/results/TMP` |
| `SLURM_TMPDIR` | SLURM job temporary directory (should match `TMPDIR`) | `$TMPDIR` |

**Why set `TMPDIR`?**

Setting `TMPDIR` explicitly is **critical** on many HPC systems:

1. **Memory-backed `/tmp` (tmpfs)**: Some HPC nodes have no local disk storage. The default `/tmp` is mounted as `tmpfs`, which stores files **in RAM**. Large temporary files from tools like STAR, RepeatMasker, or Augustus can quickly exhaust memory and crash your jobs with cryptic "out of memory", "no space left on device", or even **segmentation fault** errors—since the OS may kill processes or corrupt memory when tmpfs fills up.

2. **Quota limits**: Shared `/tmp` partitions often have strict per-user quotas (e.g., 1–10 GB). Genome annotation tools easily exceed this.

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

### Alternative: Local Execution (no SLURM)

If you don't have access to a SLURM cluster, use `annotate_local.sh` instead:

```bash
export SYLVAN_CONFIG="toydata/config/config_annotate_local.yml"

# Dry run
snakemake -n --snakefile bin/Snakefile_annotate

# Run locally (16 cores)
./bin/annotate_local.sh
```

This uses Snakemake's `--cores 16` for local parallelism instead of `--cluster`. The local config (`config_annotate_local.yml`) has per-rule thread counts scaled for a single machine (16 cores, 62 GB RAM). Adjust these values to match your hardware.

See [Local Execution without SLURM](#local-execution-without-slurm) for a full walkthrough.

**Output locations**
- All intermediate/final results are written under the repo root `results/`.
- RepeatMasker/RepeatModeler run inside `results/GETA/RepeatMasker/...`, so `.RepeatMaskerCache` and `RM_*` temp folders also stay there.
- EVM commands and outputs live in `results/EVM/`; no `EVM -> results/EVM` symlink is needed.
- For the filter pipeline, set `RexDB` to a RepeatExplorer protein DB (e.g., Viridiplantae_v4.0.fasta from https://github.com/repeatexplorer/rexdb). You can download directly via:
  `wget -O toydata/misc/Viridiplantae_v4.0.fasta https://raw.githubusercontent.com/repeatexplorer/rexdb/refs/heads/main/Viridiplantae_v4.0.fasta`

**Expected runtime (toy data):** 4–8 hours wall-clock on the [test environment](#test-environment). See [Runtime Benchmark](#runtime-benchmark) for per-step details.

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

The filter typically converges within 3–5 iterations. Each iteration adds high-confidence predictions back into the training set and retrains. The process stops when no new observations exceed the `--recycle` threshold or `--max-iter` is reached.

### Important: `chrom_regex`

The `chrom_regex` field in `config_filter.yml` is required for proper chromosome identification. It must match your genome's chromosome naming convention:

```yaml
# Common patterns:
chrom_regex: (^Chr)|(^chr)|(^LG)|(^Ch)|(^\d)
```

### Output Files

| File | Description |
|------|-------------|
| `results/FILTER/filtered.gff3` | Kept gene models |
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

## Step 5c: Feature Importance Analysis

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
  to count how many BUSCOs remain in the "keep" set during each iteration—no new
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

## Step 5d: Benchmark (BUSCO + OMArk)

Sylvan includes a benchmarking workflow that evaluates multiple GFF3 annotations
side-by-side using **BUSCO** (protein mode) and optionally **OMArk** (proteome
quality). This helps compare the quality of individual evidence sources (Augustus,
Helixer, LiftOff, GETA, EVM) against the final Sylvan output (pre-filter and
post-filter).

### Configuration

The benchmark is configured in `config_filter.yml` (or its local equivalent):

```yaml
Benchmark:
  omark_db: /usr/local/src/omark_db/LUCA.h5   # OMAmer database. Leave empty to skip OMArk.
  gff3_files:
    tair10_chr4: toydata/backup/ath4.gff3      # TAIR10 reference (ground truth)
    augustus: results/GETA/Augustus/augustus.gff3
    helixer: results/AB_INITIO/Helixer/helixer.gff3
    liftoff: results/LIFTOVER/LiftOff/liftoff.gff3
    geta: results/GETA/geta.geneModels.gff3
    geta_best: results/GETA/geta.bestGeneModels.gff3
    miniprot: results/PROTEIN/merged_rmdup_proteins.fasta.miniprot.clean.gff3
    evm: results/EVM/EVM.all.gff3
    sylvan_prefilter: results/PREFILTER/Sylvan.gff3
    sylvan_filtered: results/FILTER/filtered.gff3
```

Add or remove labels as needed. Each label maps to a GFF3 file that will be
evaluated.

### Running the benchmark

**SLURM cluster:**
```bash
./bin/benchmark.sh
```

**Local execution (no SLURM):**
```bash
./bin/benchmark_local.sh
```

**Combined pipeline (annotate + filter + benchmark):**
```bash
./bin/run_local.sh          # runs all three phases
```

You can skip individual phases with environment variables:
```bash
SYLVAN_SKIP_ANNOTATE=1 SYLVAN_SKIP_FILTER=1 ./bin/run_local.sh   # benchmark only
```

### What the benchmark does

1. **gff2pep** — Extracts protein sequences from each GFF3 using `gffread`,
   cleans headers, and strips internal stop codons.
2. **BUSCO** — Runs BUSCO in protein mode against the configured lineage
   (e.g., `eudicots_odb10`).
3. **OMArk** (optional) — Runs `omamer search` followed by `omark` for
   proteome completeness and contamination assessment.
4. **Summary** — Produces `results/BENCHMARK/benchmark_summary.tsv` with
   gene counts, BUSCO scores (C/S/D/F/M), and OMArk metrics for all labels.

### Output

```
results/BENCHMARK/
  benchmark_summary.tsv          # Main comparison table
  augustus.pep / .busco/          # Per-label protein and BUSCO results
  helixer.pep / .busco/
  ...
  sylvan_filtered.pep / .busco/
```

### OMArk setup

OMArk requires an OMAmer database file (e.g., `LUCA.h5`, ~6 GB). The pre-built
Singularity image includes it at `/usr/local/src/omark_db/LUCA.h5`. If building
your own container or running externally, download from the
[OMAmer database page](https://omabrowser.org/oma/current/) and set
`Benchmark.omark_db` in the config. If left empty, OMArk steps are skipped and
only BUSCO is run.

### Toydata benchmark results

Running the benchmark on *A. thaliana* Chr4 toydata (eudicots_odb10, 16 cores, 62 GB RAM):

**Runtime:**

| Step | Wall-clock time | Notes |
|------|----------------|-------|
| gff2pep (10 labels) | < 1 s | Parallel protein extraction |
| BUSCO (10 labels) | ~5 min | 4 parallel, protein mode |
| omamer search (10 labels) | ~3 min | 4 parallel, LUCA.h5 database |
| omark (10 labels) | ~10 min | Serialized (ete3 SQLite lock) |
| **Total (BUSCO only)** | **~5 min** | |
| **Total (BUSCO + OMArk)** | **~20 min** | |

**BUSCO results:**

| Label | Genes | C% | S% | D% | F% | M% | n |
|-------|------:|---:|---:|---:|---:|---:|--:|
| **TAIR10 (reference)** | **4,128** | **16.6** | **12.7** | **4.0** | **0.3** | **83.1** | **2,326** |
| augustus | 4,444 | 16.6 | 15.3 | 1.3 | 0.3 | 83.1 | 2,326 |
| helixer | 4,157 | 16.5 | 15.8 | 0.7 | 0.5 | 83.0 | 2,326 |
| liftoff | 11,149 | 16.7 | 9.8 | 6.9 | 0.3 | 83.0 | 2,326 |
| geta | 4,186 | 16.6 | 15.3 | 1.3 | 0.3 | 83.1 | 2,326 |
| geta_best | 4,186 | 16.6 | 15.9 | 0.7 | 0.3 | 83.1 | 2,326 |
| miniprot | 0 | 16.8 | 0.3 | 16.5 | 0.4 | 82.8 | 2,326 |
| evm | 10,625 | 7.3 | 6.9 | 0.3 | 0.4 | 92.3 | 2,326 |
| sylvan_prefilter | 5,392 | 8.0 | 7.2 | 0.8 | 0.5 | 91.6 | 2,326 |
| sylvan_filtered | 2,256 | 8.0 | 7.2 | 0.8 | 0.4 | 91.6 | 2,326 |

**OMArk results:**

| Label | Consistent% | Inconsistent% | Contamination% | Unknown% |
|-------|------------:|--------------:|---------------:|---------:|
| **TAIR10 (reference)** | **95.22** | **0.54** | **0.00** | **4.24** |
| augustus | 92.10 | 1.43 | 0.00 | 6.47 |
| helixer | 95.02 | 0.48 | 0.00 | 4.50 |
| liftoff | 94.37 | 0.57 | 0.00 | 5.05 |
| geta | 94.58 | 0.83 | 0.00 | 4.59 |
| geta_best | 94.51 | 0.84 | 0.00 | 4.66 |
| miniprot | 95.77 | 1.23 | 0.00 | 3.00 |
| evm | 45.66 | 2.94 | 0.00 | 51.40 |
| sylvan_prefilter | 48.87 | 2.83 | 0.00 | 48.30 |
| sylvan_filtered | 84.91 | 1.01 | 0.00 | 14.08 |

> **Notes:**
> - **TAIR10 (reference)** is the Phytozome v9.0 annotation for Chr4 (4,128 protein-coding
>   genes). This serves as the ground truth for comparison.
> - Low BUSCO C% is expected for single-chromosome data — only ~17% of eudicot BUSCOs
>   map to Chr4. The reference TAIR10 itself shows 16.6% C, confirming this ceiling.
> - **EVM/prefilter BUSCO drop (16.6% → 7.3%)**: EVM produces many fragmented gene models
>   without proper start codons (only 36% of EVM proteins start with Met). This is a
>   toydata-specific artifact caused by gene fragmentation at chromosome boundaries during
>   EVM partitioning. On whole-genome runs, EVM scores are comparable to individual evidence
>   sources.
> - **Filter aggressiveness**: The filter reduces prefilter (5,392) to filtered (2,256) —
>   a 58% reduction. This is more aggressive than expected because many EVM-derived models
>   lack homology/expression evidence due to fragmentation, causing the RF classifier to
>   discard them. On whole-genome data with higher-quality EVM models, the retention rate
>   is typically 60-70%.
> - High "Unknown" in EVM/prefilter OMArk reflects partial gene models that lack
>   OMAmer hits. Filtering (sylvan_filtered) removes these, raising Consistent% from
>   49% to 85%.
> - 0% Contamination confirms all gene models map to *Arabidopsis thaliana*.

## Step 6: Format Output (TidyGFF)

```bash
singularity exec sylvan.sif python bin/TidyGFF.py \
  Ath4 results/FILTER/filtered.gff3 \
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
- Final outputs (`PREFILTER/Sylvan.gff3`, `FILTER/`)
- Log files (`results/logs/`)
- Configuration files, Singularity images, and pipeline scripts
- LncDC database files

> **Warning:** Only run this after you have verified both pipeline phases completed successfully. The annotation phase intermediate files cannot be regenerated without re-running the full annotation pipeline.

---

## Advanced Configuration

### EVM Weights (`evm_weights.txt`)

Controls how EvidenceModeler prioritizes evidence sources. See [README — EVM Weights](README.md#evm-weights-evm_weightstxt) for the weight table and tuning guidance.

### EVM Partitioning (`num_evm_files`)

The `num_evm_files` parameter in `config_annotate.yml` controls how many parallel EVM partitions to create. The genome is split into overlapping segments (1 Mb segments with 20 kb overlap):

- **Higher values** = more parallel SLURM jobs = faster wall-clock time, but more job scheduling overhead
- **Lower values** = fewer jobs = less cluster burden, but slower
- Default: `126` (works well for genomes up to ~500 Mb)
- For very large genomes (>1 Gb), consider increasing to 200–500

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

### Custom Helixer Models

To use custom Helixer `.h5` model files instead of the container defaults:

1. Set `helixer_model_dir` in your config to the host directory containing the model files:
   ```yaml
   helixer_model_dir: "/path/to/custom/models"
   ```

2. Bind the directory into the container via `SINGULARITY_BIND`:
   ```bash
   export SINGULARITY_BIND="/path/to/custom/models"
   ```

The pipeline will look for `{helixer_model_dir}/{helixer_model}.h5` (e.g., `/path/to/custom/models/land_plant.h5`). Leave `helixer_model_dir` empty (default) to use the container's built-in models at `/usr/local/src/Helixer/models/`.

### Customizing SLURM Submission

> **Tip:** SLURM settings can live in a separate cluster YAML (`cluster_annotate.yml` / `cluster_filter.yml`) or in the pipeline config's `__default__` section (single-file mode). Set `SYLVAN_CLUSTER_CONFIG` to point to a separate cluster file, or leave it unset to use the pipeline config.

Job submission is handled by `bin/cluster_submit.py`, which dynamically builds the `sbatch` command from the `__default__` section of your cluster config. **Account (`-A`) and partition (`-p`) are automatically skipped** when set to empty or `placeholder`.

```yaml
__default__:
  account: your-account       # Leave empty or "placeholder" if not required
  partition: your-partition    # Leave empty or "placeholder" if not required
  time: "3-00:00:00"           # Always quote time values (YAML 1.1 treats unquoted HH:MM:SS as numbers)
  memory: 8g
  extra_args: ""              # Extra sbatch flags (e.g., "--export=ALL")
```

> **Warning:** Always quote YAML time values (e.g., `time: "3-00:00:00"`). Unquoted values like `72:00:00` are silently parsed as sexagesimal integers by YAML 1.1, resulting in incorrect SLURM walltimes.

For per-rule customization (e.g., GPU for Helixer), use `extra_args` in the rule section:
```yaml
helixer:
  extra_args: "--gres=gpu:1"
```

See [README — Customizing SLURM Submission](README.md#customizing-slurm-submission) for more details.

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

For the full troubleshooting table, see [README — Common Issues](README.md#common-issues).

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

**Annotation Phase (PREFILTER/Sylvan.gff3):**

| Metric | Count |
|--------|-------|
| Total genes | 5,720 |
| Total mRNA | 5,800 |

**Filter Phase (filtered.gff3):**

| Metric | Count |
|--------|-------|
| Genes kept | 3,756 |
| mRNA kept | 3,834 |
| Genes discarded | 1,964 |
| mRNA discarded | 1,966 |

**Output files:**
- `results/FILTER/filtered.gff3` - Kept gene models
- `results/FILTER/discard.gff3` - Discarded gene models
- `results/FILTER/data.tsv` - Feature matrix used by random forest (input to feature importance analysis)
- `results/FILTER/keep_data.tsv` - Evidence data for kept genes
- `results/FILTER/discard_data.tsv` - Evidence data for discarded genes
- `results/FILTER/{prefix}.cdna` - Extracted transcript sequences
- `results/FILTER/{prefix}.pep` - Extracted peptide sequences
- `results/FILTER/lncrna_predict.csv` - lncDC long non-coding RNA predictions
- `results/PREFILTER/Sylvan.gff3.map` - ID mapping between original and new IDs

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

---

## Utility Scripts

### Generate Cluster Config

`bin/generate_cluster_from_config.py` extracts the SLURM-relevant sections (`__default__` and per-rule overrides) from `config_annotate.yml` into a standalone cluster-only YAML. This is optional — `config_annotate.yml` already serves as both pipeline config and `--cluster-config`, so you only need this script if you want a minimal, separate cluster file.

```bash
# Extract standalone cluster YAML from production config
# Default: auto-detects partition max time from sinfo and sets time = max - 1 day
python bin/generate_cluster_from_config.py \
  --config config/config_annotate.yml \
  --out config/cluster_annotate.yml \
  --account your-account --partition your-partition

# Extract from toydata config (explicit time override)
python bin/generate_cluster_from_config.py \
  --config toydata/config/config_annotate.yml \
  --out toydata/config/cluster_annotate.yml \
  --account your-account --partition your-partition \
  --time "5-00:00:00"
```

**Auto time detection** (default): The script queries `sinfo` for the partition's maximum walltime and sets `time` to `max - 1 day`. If `sinfo` is unavailable or the partition is not found, it falls back to `9-00:00:00`. Override with `--time "D-HH:MM:SS"` to set a specific walltime.

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
| `cluster_submit.py` | SLURM job submission wrapper — skips `-A`/`-p` when account/partition are empty |
| `MonitorFilter.py` | Visualize filter iteration progress (matplotlib) |

---

## Runtime Benchmark

The following benchmarks were measured on the toy dataset (*A. thaliana* chromosome 4, 18.6 Mb) using the test environment described below.

### Test Environment

| Specification | Value |
|---------------|-------|
| Nodes | 4 |
| Total CPUs | 256 |
| CPU | Intel Xeon E5-2683 v4 @ 2.10GHz |
| Cores per node | 64 (2 sockets x 16 cores x 2 threads) |
| Memory per node | 256 GB |
| Storage | GPFS |
| GPU | NVIDIA V100 (Helixer only) |

### Key Observations

- **Most time-consuming steps**: `aggregate_CombineGeneModels` (~50,000 s), `geneRegion2Genewise` (~1,000 s), and `Sam2Transfrag` (~100–200 s) are the bottlenecks
- **Fast steps**: Most preprocessing and formatting rules complete in under 10 seconds
- **Parallelizable rules**: Rules like `geneRegion2Genewise`, `gmapExon`, and `STAR_paired` run as multiple parallel jobs, significantly reducing wall-clock time
- **GPU-accelerated**: `helixer` benefits from GPU acceleration when available

### Runtime Variability

Actual runtime will vary significantly depending on your hardware and cluster configuration:

| Factor | Impact |
|--------|--------|
| **CPU speed** | Faster clock speeds reduce single-threaded bottlenecks |
| **Available nodes** | More nodes = more parallel jobs = faster wall-clock time |
| **Memory per node** | Insufficient memory causes job failures or swapping |
| **Storage I/O** | GPFS/Lustre faster than NFS; SSD faster than HDD |
| **Queue wait time** | Busy clusters add significant delays between jobs |
| **GPU availability** | Helixer runs ~10x faster with GPU acceleration |

### Phase 1 — Annotate

| Step | Wall-clock time | Notes |
|------|----------------|-------|
| EDTA (repeat library) | ~1.5 h | Pre-pipeline; LINE detection is the bottleneck |
| RepeatMasking | 15–30 min | |
| RNA-seq alignment (STAR) | 30–60 min | Parallelized across samples |
| Transcript assembly (StringTie/PsiCLASS) | 20–40 min | |
| Protein homology (Miniprot/GeneWise) | 30–60 min | Parallelized across genome regions |
| **Helixer (GPU)** | **~5 min** | **Single NVIDIA V100; see GPU note below** |
| Augustus training + prediction | 1–2 h | Longest single-threaded bottleneck |
| Gene model combination (CombineDuck) | 30–60 min | ~50,000 s total CPU across parallel jobs |
| EvidenceModeler | 30–60 min | |
| **Phase 1 total** | **4–8 h** | Wall-clock with 4 nodes / 256 CPUs |

### Phase 2 — Filter

| Step | Wall-clock time | Notes |
|------|----------------|-------|
| CDS/peptide extraction | < 1 min | |
| PfamScan | 5–10 min | |
| RSEM (expression quantification) | 10–20 min | |
| BLASTp + RexDB | 10–20 min | Parallelized across 20 splits |
| Ab initio coverage (Helixer/Augustus/RM) | < 5 min | |
| lncDC | < 5 min | |
| BUSCO | 5–10 min | |
| Semi-supervised random forest | < 5 min | 3–5 iterations to convergence |
| **Phase 2 total** | **~1 h** | |

### End-to-End Summary

| Component | Time |
|-----------|------|
| EDTA (pre-pipeline) | ~1.5 h |
| Phase 1 — Annotate | 4–8 h |
| Phase 2 — Filter | ~1 h |
| **Total** | **~7–11 h** |

With the test environment above (4 nodes, 256 CPUs, 256 GB/node), the toy dataset completes in **4–8 hours** wall-clock time for Phase 1 alone. On smaller clusters or shared resources, expect longer runtimes.

### Helixer GPU vs CPU

Helixer is the only pipeline step that benefits from GPU acceleration. On the toy dataset (18.6 Mb):

| Hardware | Helixer runtime | Notes |
|----------|----------------|-------|
| NVIDIA V100 (1 GPU) | ~5 min | Default in pipeline; `--gres=gpu:1` in SLURM config |
| CPU-only (16 threads) | ~50 min | ~10x slower; set via Helixer `--no-gpu` flag |

Helixer GPU time represents < 2% of total Sylvan wall-clock time, so GPU availability is helpful but not a bottleneck. The pipeline runs end-to-end without a GPU if needed.

### Scalability for EBP-Scale Genomes

The Earth BioGenome Project (EBP) targets thousands of species with genomes ranging from hundreds of Mb to several Gb. Key scaling considerations:

| Factor | Scaling behavior |
|--------|-----------------|
| **Genome size** | Most steps scale linearly (repeat masking, alignment, EVM). Augustus training and CombineDuck scale super-linearly due to increased model complexity. |
| **RNA-seq volume** | STAR/HiSat2 alignment scales linearly with read count. Parallelized across samples. |
| **Protein DB size** | Miniprot/GeneWise scale linearly. BLASTp is parallelized across 20 splits. |
| **Helixer (GPU)** | Scales linearly with genome size via `--subsequence-length`. A 1 Gb genome takes ~30–60 min on a V100. |
| **Filter phase** | Scales with gene count, not genome size. Typically < 2 h even for large annotations (50k+ genes). |

**Estimated runtimes by genome size** (4 nodes / 256 CPUs / 1 GPU):

| Genome size | Example organism | Estimated wall-clock |
|-------------|-----------------|---------------------|
| 18.6 Mb (toy) | *A. thaliana* chr4 | 7–11 h |
| ~135 Mb | *A. thaliana* (full) | 1–2 days |
| ~750 Mb | Rice, tomato | 3–5 days |
| ~2.5 Gb | Maize, wheat subgenome | 1–2 weeks |

These estimates assume continuous job scheduling. In practice, SLURM queue wait times can dominate wall-clock time on shared clusters. Snakemake's job-level parallelism means that adding nodes provides near-linear speedup for the parallelizable steps (RNA-seq alignment, GeneWise, BLASTp, EVM).

---

## Local Execution (without SLURM)

This section describes how to run the full Sylvan pipeline on a single Linux machine without SLURM.

### Local Test Environment

| Specification | Value |
|---------------|-------|
| CPUs | 16 cores |
| Memory | 62 GB RAM |
| GPU | NVIDIA GPU (CUDA 12.6) |
| Storage | Local SSD |
| OS | Ubuntu 22.04 (Linux 6.8) |
| Singularity | 3.x (writable sandbox) |

### Step-by-Step

**1. Create a local config**

Copy the provided local config template:

```bash
cp toydata/config/config_annotate_local.yml my_local_config.yml
```

Key differences from the SLURM config:
- `account` and `partition` under `__default__` are empty (ignored)
- Per-rule `ncpus`/`threads` are capped for a single machine (e.g., max 12 for heavy rules)

**2. Run the pipeline**

```bash
export SYLVAN_CONFIG="toydata/config/config_annotate_local.yml"

# Dry run first
snakemake -n --snakefile bin/Snakefile_annotate

# Run
./bin/annotate_local.sh
```

**3. Monitor progress**

```bash
# Snakemake log
tail -f .snakemake/log/*.snakemake.log

# Per-rule logs
ls -lt results/logs/*.err | head -10
```

### Local Execution: Issues Encountered and Fixes

During local testing on the toy dataset, 17 issues were identified and fixed. These are documented in [`error.md`](error.md). Key categories:

**Container execution issues:**
- `run:` blocks in Snakemake execute on the HOST, not inside the container. Fixed with a `run_in_container()` helper that wraps commands with `singularity exec`.
- Affected rules: `mergeTransfrag`, `geneRegion2Genewise`, `geneModels2AugusutsTrainingInput`, `BGM2AT`, `augustusWithHints`

**Shell compatibility (dash vs bash):**
- The container's `/bin/sh` is dash (Ubuntu-based), not bash. Perl `system()` calls use `/bin/sh`.
- `&>` in dash means "background + redirect", not "redirect stderr". Fix: use `> file 2>&1`.

**Tool-specific fixes:**
- RepeatModeler 2.0.5: `-pa` deprecated → use `-threads`; output path changed to `RM_*/consensi.fa`
- TransDecoder 5.7.1: `.gff3` output may not be retained → fallback recovery from internal checkpoints
- Augustus config: `cp -rf config/ target/` creates nested directory → use `cp -rf config/* target/`
- AGAT: `-gff` parsed as `-g ff` (Getopt::Long bundling) → use `--gff`
- `hmmscan` in `filter` env, not `genepred`

**GPU/CUDA:**
- Container v4+ bundles TensorFlow with pip CUDA packages (`nvidia-cuda-runtime-cu12`, `nvidia-cudnn-cu12`), eliminating host CUDA dependency. Only the NVIDIA driver (>= 525.60.13) is needed on the host. The `--nv` flag is safe on CPU-only nodes (Singularity silently skips it).

### Local Runtime Benchmark (Toy Data)

Measured on the local test environment above (16 cores, 62 GB RAM, single machine) with
Sylvan container v3 (CUDA 12.2):

| Phase | Wall-clock time | Steps | Notes |
|-------|----------------|-------|-------|
| Phase 1 — Annotate | ~8 hours | 228 | 16 cores, Helixer + Augustus + EVM |
| Phase 2 — Filter | ~22 min | 94 | RF filtering, BUSCO, PfamScan |
| Phase 3 — Benchmark (BUSCO only) | ~5 min | 28 | BUSCO protein mode, 10 GFF3 files |
| Phase 3 — Benchmark (BUSCO + OMArk) | ~20 min | 46 | + omamer search + omark (serialized) |
| **Total (BUSCO only)** | **~8.5 hours** | **350** | |
| **Total (BUSCO + OMArk)** | **~8.7 hours** | **368** | |

> Helixer now produces real output with the CUDA 12.2 container (previously
> produced 0 genes with CUDA 11.2). See [Step 5d](#step-5d-benchmark-busco--omark)
> for full BUSCO and OMArk results.

**Key outputs (local run with v3 container):**
- `results/PREFILTER/Sylvan.gff3` — 7.8 MB, 5,392 gene models
- `results/FILTER/filtered.gff3` — 3.8 MB, 2,256 gene models
- `results/BENCHMARK/benchmark_summary.tsv` — comparison table

---

## Getting Help

- Issues: https://github.com/plantgenomicslab/Sylvan/issues
- See also: [README.md](README.md) for configuration reference
