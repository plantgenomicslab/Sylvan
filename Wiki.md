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
