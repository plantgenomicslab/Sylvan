# Sylvan Tutorial: Running with Toy Data

This tutorial walks you through running Sylvan on the included toy dataset (*A. thaliana* chromosome 4).

## Table of Contents

- [Prerequisites](#prerequisites)
- [Step 1: Environment Setup](#step-1-environment-setup)
- [Step 2: Download Containers](#step-2-download-containers)
- [Step 3: Prepare Repeat Library (EDTA)](#step-3-prepare-repeat-library-edta)
- [Step 4: Run Annotation Pipeline](#step-4-run-annotation-pipeline)
- [Step 5: Run Filter Pipeline](#step-5-run-filter-pipeline)
- [Step 6: Format Output (TidyGFF)](#step-6-format-output-tidygff)
- [Monitoring and Debugging](#monitoring-and-debugging)
- [Toy Data Details](#toy-data-details)

---

## Prerequisites

- Linux system with SLURM
- Singularity 3.x+
- Conda/Mamba
- Git LFS

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

## Step 4: Run Annotation Pipeline

### Environment Variables

The pipeline uses several environment variables for configuration:

| Variable | Description | Example |
|----------|-------------|---------|
| `SYLVAN_CONFIG` | Path to pipeline config file | `toydata/config/config_annotate.yml` |
| `SYLVAN_CLUSTER_CONFIG` | Path to SLURM cluster config (optional, auto-derived from `SYLVAN_CONFIG`) | `toydata/config/cluster_annotate.yml` |
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

After annotation completes:

```bash
sbatch -A [account] -p [partition] -c 1 --mem=4g \
  -J filter -o filter.out -e filter.err \
  --wrap="./bin/filter.sh"
```

**Output:** `results/FILTER/filter.gff3`

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
| Out of memory | Increase `memory` in `cluster_annotate.yml` |
| LFS files are pointers | Run `git lfs pull` |
| Singularity bind error | Use paths within working directory |
| Augustus training fails | Need minimum 500 training genes |

---

## Toy Data Details

### Overview

The toy dataset contains *Arabidopsis thaliana* chromosome 4 split into 3 segments (~18.6 Mb total).

### Directory Structure

```
toydata/
├── config/                      # Configuration files
│   ├── config_annotate.yml      # Pipeline config (pre-configured)
│   ├── cluster_annotate.yml     # SLURM resource config
│   └── evm_weights.txt          # EVM weights
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

> **Note:** The original chromosome 4 paper reported 3,744 protein-coding genes. TAIR10 increased this count to 4,124 as annotation methods improved.
>
> Source: [Phoenix Bioinformatics / TAIR](https://www.arabidopsis.org/)

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

The total wall-clock time depends heavily on cluster availability and parallelization. With sufficient resources, the toy dataset completes in 4-8 hours.

---

## Getting Help

- Issues: https://github.com/plantgenomicslab/Sylvan/issues
- See also: [README.md](README.md) for configuration reference
