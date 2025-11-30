# Sylvan Tutorial: Running with Toy Data

![Sylvan Workflow](https://github.com/plantgenomicslab/Sylvan/blob/main/docs/images/sylvan_workflow.jpg?raw=true)

This tutorial guides you through setting up the environment, downloading the Singularity image, and running the EDTA step using toy data.

## 1. Set up Conda Environment

Create and activate the `sylvan` environment with the necessary dependencies. We prioritize `conda-forge` to ensure dependency resolution.

```bash
# Create the environment
conda create -n sylvan -c conda-forge -c bioconda \
  python=3.11 \
  snakemake=7 -y

# Activate the environment
conda activate sylvan
```

## 2. Download Singularity Image

Pull the Sylvan Singularity image from [Sylabs Cloud Library](https://cloud.sylabs.io/library/wyim/sylvan/sylvan).

```bash
singularity pull --arch amd64 sylvan.sif library://wyim/sylvan/sylvan:latest
```

## Git LFS (Large File Storage)

The toy data files are stored using [Git LFS](https://git-lfs.github.com/) to keep the repository lightweight. You must have Git LFS installed to properly clone the repository with all toy data files.

### Installing Git LFS

```bash
# Ubuntu/Debian
sudo apt-get install git-lfs

# macOS (Homebrew)
brew install git-lfs

# CentOS/RHEL
sudo yum install git-lfs

# Initialize Git LFS (run once after installation)
git lfs install
```

### Cloning with LFS

When you clone the repository, Git LFS files are automatically downloaded if Git LFS is installed:

```bash
git clone https://github.com/plantgenomicslab/Sylvan.git
```

If you already cloned the repository without Git LFS, you can fetch the LFS files:

```bash
cd Sylvan
git lfs pull
```

### Verifying LFS Files

To check if LFS files are properly downloaded (not just pointers):

```bash
# List LFS-tracked files
git lfs ls-files

# Check file sizes (should be > 200 bytes if downloaded)
ls -la toydata/
```

> **Note:** If toydata files appear as small text files (~130 bytes) containing `version https://git-lfs.github.com/spec/v1`, Git LFS is not properly configured. Run `git lfs pull` to download the actual files.

## Toy Data Overview

The toy dataset was constructed from *Arabidopsis thaliana* chromosome 4, split into 3 equal parts.

### Toydata Directory Structure

```
toydata/
├── config/                      # Configuration files
│   ├── config_annotate.yml      # Annotation phase config (pre-configured)
│   ├── config_filter.yml        # Filter phase config
│   └── evm_weights.txt          # EVM weights
├── genome_input/                # Input genome
│   └── genome.fasta.gz          # A. thaliana Chr4 (3 parts, ~18.6 Mb)
├── RNASeq/                      # RNA-seq data (paired-end)
│   ├── sub_SRR1019221_1.fastq.gz
│   ├── sub_SRR1019221_2.fastq.gz
│   └── ...                      # Additional paired samples
├── neighbor_genome/             # Neighbor species genomes for Liftoff
│   ├── aly4.fasta               # A. lyrata Chr4
│   ├── ath4.fasta               # A. thaliana Chr4
│   ├── chi4.fasta               # C. hirsuta Chr4
│   └── cru4.fasta               # C. rubella Chr4
├── neighbor_gff3/               # Neighbor species annotations for Liftoff
│   ├── aly4.gff3
│   ├── ath4.gff3
│   ├── chi4.gff3
│   └── cru4.gff3
├── cds_aa/                      # CDS and protein sequences
│   ├── neighbor.cds             # Combined CDS for EDTA
│   └── neighbor.aa              # Protein sequences for homology
└── EDTA/                        # Pre-computed EDTA output
    ├── genome.fasta.mod.EDTA.TElib.fa    # Repeat library
    ├── genome.fasta.mod.EDTA.TEanno.gff3 # TE annotations
    └── neighbor.cds             # CDS used for EDTA
```

### How the Toy Data Was Created

**1. Genome (genome.fasta)**

*Arabidopsis thaliana* chromosome 4 (~18.6 Mb) was divided into three equal segments:
- Chr4_1 (6,195,060 bp)
- Chr4_2 (6,195,060 bp)
- Chr4_3 (6,195,018 bp)

```bash
# Split chromosome 4 into 3 equal parts
seqkit split2 -p 3 Chr4.fasta
```

**2. Neighbor CDS (neighbor.cds)**

To provide CDS evidence for EDTA repeat annotation, syntenic regions were identified between *A. thaliana* Chr4 and closely related Brassicaceae species using [MCscan (jcvi)](https://github.com/tanghaibao/jcvi) ([Tang et al., 2024, iMeta](https://doi.org/10.1002/imt2.211)).

Related species used for synteny analysis:
- *Arabidopsis lyrata*
- *Capsella rubella*
- *Cardamine hirsuta*

```bash
# Run MCscan to identify syntenic blocks
python -m jcvi.compara.catalog ortholog ath aly --no_strip_names
python -m jcvi.compara.catalog ortholog ath cru --no_strip_names
python -m jcvi.compara.catalog ortholog ath chi --no_strip_names

# Extract CDS from syntenic regions
python -m jcvi.compara.synteny mcscan ath.bed ath.aly.lifted.anchors --iter=1 -o ath.aly.i1.blocks
```

CDS sequences from genes within syntenic blocks were extracted and combined into `neighbor.cds`.

**3. Repeat Library**

EDTA 2.2.0 was run on the toy genome with the neighbor CDS to generate the repeat library.

Download the EDTA container (released 2024-04):
```bash
export SINGULARITY_CACHEDIR=$PWD
singularity pull EDTA.sif docker://quay.io/biocontainers/edta:2.2.0--hdfd78af_1
```

Run EDTA:
```bash
singularity exec --cleanenv --env PYTHONNOUSERSITE=1 EDTA.sif EDTA.pl \
    --genome genome.fasta \
    --cds neighbor.cds \
    --anno 1 \
    --threads 16
```

Output: `genome.fasta.mod.EDTA.TElib.fa`

**4. RNA-seq Data**

*A. thaliana* RNA-seq samples were downloaded from NCBI SRA, quality-filtered with fastp, and mapped to Chr4 using STAR. Only reads mapping to Chr4 were extracted to create a lightweight test dataset:

```bash
# Download from SRA
prefetch SRR1019221
fasterq-dump SRR1019221

# Map to Chr4 with STAR
STAR --genomeDir star_index \
     --readFilesIn SRR1019221_1.fastq.gz SRR1019221_2.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate

# Extract only mapped reads
samtools view -b -F 4 Aligned.sortedByCoord.out.bam > mapped.bam
samtools fastq -1 SRR1019221_1.fastq -2 SRR1019221_2.fastq mapped.bam
```

The following SRA samples were processed:

| SRA Accession | Tissue | Condition | Platform | Size |
|---------------|--------|-----------|----------|------|
| SRR1019221 | Leaf (14-day) | Col-0 WT | HiSeq 2000 | 4.6 Gb |
| SRR1105822 | Rosette (19-day) | Col-0 WT | HiSeq 2500 | 3.1 Gb |
| SRR1105823 | Rosette (19-day) | Col-0 WT | HiSeq 2500 | 4.5 Gb |
| SRR1106559 | Rosette (19-day) | Col-0 WT | HiSeq 2500 | 3.6 Gb |
| SRR446027 | Whole plant | Mock 1hpi (replicate A) | GA II | 2.5 Gb |
| SRR446028 | Whole plant | Mock 1hpi (replicate B) | GA II | 5.4 Gb |
| SRR446033 | Whole plant | Mock 6hpi (replicate A) | GA II | 5.3 Gb |
| SRR446034 | Whole plant | Mock 6hpi (replicate B) | GA II | 5.4 Gb |
| SRR446039 | Whole plant | Mock 12hpi (replicate A) | GA II | 2.5 Gb |
| SRR446040 | Whole plant | Mock 12hpi (replicate B) | GA II | 6.3 Gb |
| SRR764885 | Leaf (4-week) | Col-0 WT | HiSeq 2000 | 4.6 Gb |
| SRR934391 | Whole plant | Col-0 WT | GA IIx | 4.0 Gb |

### Sequence Statistics

| ID | Length | GC (%) | ΔGC |
| :--- | :--- | :--- | :--- |
| Chr4_1 | 6195060 | 36.66 | -0.30 |
| Chr4_2 | 6195060 | 35.24 | -0.52 |
| Chr4_3 | 6195018 | 36.69 | 0.13 |

### File Summary

| file | format | type | num_seqs | sum_len | min_len | avg_len | max_len | Q1 | Q2 | Q3 | sum_gap | N50 | Q20(%) | Q30(%) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| genome.fasta | FASTA | DNA | 3 | 18585138 | 6195018 | 6195046.0 | 6195060 | 6195018 | 6195060 | 6195060 | 0 | 6195060 | 0.00 | 0.00 |

### Toydata Directory Structure

```
toydata/
├── genome.fasta              # Toy genome (A. thaliana Chr4 split into 3 parts)
├── neighbor.cds              # Combined CDS from all neighbor species
│
├── aly4.fasta                # A. lyrata Chr4 genome
├── aly4.gff3                 # A. lyrata Chr4 annotation
├── aly4.cds                  # A. lyrata Chr4 CDS sequences
├── aly4.ids                  # A. lyrata gene IDs
├── aly4.fasta.fai            # A. lyrata FASTA index
│
├── cru4.fasta                # C. rubella Chr4 genome
├── cru4.gff3                 # C. rubella Chr4 annotation
├── cru4.cds                  # C. rubella Chr4 CDS sequences
├── cru4.ids                  # C. rubella gene IDs
├── cru4.fasta.fai            # C. rubella FASTA index
│
├── chi4.fasta                # C. hirsuta Chr4 genome
├── chi4.gff3                 # C. hirsuta Chr4 annotation
├── chi4.cds                  # C. hirsuta Chr4 CDS sequences
├── chi4.ids                  # C. hirsuta gene IDs
├── chi4.fasta.fai            # C. hirsuta FASTA index
│
├── ath4.fasta                # A. thaliana Chr4 genome
├── ath4.gff3                 # A. thaliana Chr4 annotation
├── ath4.cds                  # A. thaliana Chr4 CDS sequences
├── ath4.ids                  # A. thaliana gene IDs
├── ath4.fasta.fai            # A. thaliana FASTA index
│
└── [STAR index files]        # Pre-built STAR index
    ├── Genome
    ├── SA
    ├── SAindex
    ├── chrLength.txt
    ├── chrName.txt
    ├── chrNameLength.txt
    ├── chrStart.txt
    └── genomeParameters.txt
```

**Species abbreviations:**
| Code | Species | Common Name |
|------|---------|-------------|
| aly4 | *Arabidopsis lyrata* | Lyrate rockcress |
| cru4 | *Capsella rubella* | Pink shepherd's purse |
| chi4 | *Cardamine hirsuta* | Hairy bittercress |
| ath4 | *Arabidopsis thaliana* | Thale cress |

## Test Environment

The toy data was tested on a high-performance computing cluster with the following specifications:

### Cluster Overview

| Partition | Nodes | CPU Cores | Memory | GPU |
|-----------|-------|-----------|--------|-----|
| CPU | 108 | 3,456 | 24.8 TiB | - |
| GPU | - | 352 | 2.75 TiB | 44x NVIDIA Tesla P100 |
| Visualization | - | 48 | 1.1 TiB | 3x NVIDIA Tesla V100 |

### Node Specifications

| Specification | Value |
|---------------|-------|
| Architecture | x86_64 |
| CPU | Intel Xeon E5-2683 v4 @ 2.10GHz |
| Sockets | 2 |
| Cores per socket | 16 |
| Threads per core | 2 |
| Total CPUs | 64 |
| L3 Cache | 40 MB |
| CPU Max MHz | 3000 |

## 3. Download EDTA Container

For EDTA, we recommend using the official EDTA container to avoid TensorFlow compatibility issues.

```bash
export SINGULARITY_CACHEDIR=$PWD
singularity pull EDTA.sif docker://quay.io/biocontainers/edta:2.2.0--hdfd78af_1
```

## 4. Run EDTA on Toy Data

Submit a SLURM job to run EDTA on the toy dataset. Ensure your `toydata` directory contains the necessary `genome.fasta` and `neighbor.cds` files.

```bash
sbatch -c 16 --mem=68g --wrap="/usr/bin/time -v singularity exec --cleanenv --env PYTHONNOUSERSITE=1 EDTA.sif EDTA.pl --genome toydata/genome.fasta --cds toydata/neighbor.cds --anno 1 --threads 16 --force 1"
```

**Supported file extensions for neighbor files:**
- **Genome assemblies:** `.fa` or `.fasta` (NOT `.fas`)
- **Annotations:** `.gff` or `.gff3`

> **Note:** The `--force 1` parameter is required for this toy dataset because:
> 1. **Subset genome**: The toy data contains only a small portion of the *A. thaliana* genome (~18.6 Mb), which may lack intact TIR elements
> 2. **TIR-Learner3.0 compatibility**: TensorFlow/Keras version incompatibility (`AttributeError: as_proto`) may cause TIR detection to fail
>
> This option allows EDTA to continue even when certain TE types are not detected, producing results for LTR, SINE, LINE, and Helitron elements.

### Benchmark Results

| Stage | Duration | Notes |
|-------|----------|-------|
| LTR detection | ~3 min | LTR_retriever + LTR_FINDER |
| SINE detection | ~6 min | AnnoSINE |
| LINE detection | ~45 min | May produce 0 bp for small genomes |
| TIR detection | ~8 min | TIR-Learner3.0 (may fail, use --force 1) |
| Helitron detection | ~10 min | HelitronScanner |
| Filtering & annotation | ~8 min | TEsorter + RepeatMasker |
| **Total** | **~1.5 hours** | 18.6 Mb genome, 16 threads |

### Expected Output

```
genome.fasta.mod.EDTA.TElib.fa      # Final TE library for RepeatMasker
genome.fasta.mod.EDTA.TEanno.gff3   # Whole-genome TE annotation
genome.fasta.mod.EDTA.TEanno.sum    # TE annotation summary
genome.fasta.mod.MAKER.masked       # Masked genome for gene annotation
```

### TE Annotation Results (by Sequence)

| Sequence | TE Count | Masked (bp) | Masked (%) |
|----------|----------|-------------|------------|
| Chr4_1 | 1,831 | 1,416,349 | 22.86% |
| Chr4_2 | 283 | 281,871 | 4.55% |
| Chr4_3 | 73 | 108,260 | 1.75% |
| **Total** | **2,187** | **1,806,480** | **9.72%** |

## 5. Run Sylvan Annotation Pipeline

After completing EDTA, you can run the main Sylvan annotation pipeline. The pipeline uses Snakemake for workflow management.

### Configuration

The toy data includes a pre-configured configuration file at `toydata/config/config_annotate.yml`. For your own data, copy and modify `config/config_annotate.yml`.

Key configuration parameters:

| Parameter | Description | Toy Data Value |
|-----------|-------------|----------------|
| `prefix` | Output file prefix | `ath4_toydata` |
| `genome` | Genome FASTA (.fa, .fasta, .fa.gz, .fasta.gz) | `toydata/genome_input/genome.fasta.gz` |
| `rna_seq` | RNA-seq data directory | `toydata/RNASeq` |
| `proteins` | Protein FASTA for homology | `toydata/cds_aa/neighbor.aa` |
| `singularity` | Singularity image path | `sylvan.sif` |
| `neighbor_gff` | Neighbor species GFF3 directory | `toydata/neighbor_gff3` |
| `neighbor_fasta` | Neighbor species genome directory | `toydata/neighbor_genome` |
| `RM_lib` | Repeat nucleotide FASTA file | `toydata/EDTA/genome.fasta.mod.EDTA.TElib.fa` |
| `evm_weights` | EVM weights file | `toydata/config/evm_weights.txt` |

### Running the Pipeline

#### Dry Run (Recommended First)

Test your configuration without executing any jobs:

```bash
snakemake -s bin/Snakefile_annotate \
    --configfile toydata/config/config_annotate.yml \
    --dry-run
```

#### Local Execution

Run on a local machine (for testing only):

```bash
snakemake -s bin/Snakefile_annotate \
    --configfile toydata/config/config_annotate.yml \
    --cores 8 \
    --use-singularity
```

#### SLURM Cluster Execution

For HPC clusters with SLURM, use the snakemake executor plugin:

```bash
snakemake -s bin/Snakefile_annotate \
    --configfile toydata/config/config_annotate.yml \
    --workflow-profile profiles/slurm \
    --jobs 100 \
    --use-singularity
```

### Output Directory Structure

All pipeline outputs are organized under the `results/` directory:

```
results/
├── AB_INITIO/          # Ab initio gene predictions (Helixer)
├── EVM/                # EVM consensus gene models
├── FILTER/             # Filtered final output (after filter phase)
├── GETA/               # GETA pipeline outputs
│   ├── RepeatMasker/
│   │   └── .RepeatMaskerCache/  # Cache for potential rerun
│   └── ...
├── LIFTOVER/           # Liftoff annotations from neighbor species
├── PROTEIN/            # Protein alignments
├── TRANSCRIPT/         # Transcript assemblies
│   ├── evigene/        # Evigene intermediate files (for rerun)
│   │   ├── okayset/
│   │   ├── dropset/
│   │   ├── inputset/
│   │   └── intermediate_files/
│   ├── PASA/
│   │   ├── intermediate/       # PASA checkpoints (for rerun)
│   │   └── intermediate_post/  # PASA post checkpoints
│   └── ...
├── TMP/                # Temporary files
├── logs/               # SLURM job logs ({rule}_{wildcards}.out/err)
├── complete_draft.gff3 # Final combined annotation
└── EVM.all.gff3        # EVM consensus predictions
```

> **Note:** Intermediate files from evigene, PASA, and RepeatMasker are preserved in `results/` for potential reruns.

### Monitoring Progress

View running jobs:
```bash
squeue -u $USER
```

Check Snakemake log:
```bash
tail -f .snakemake/log/*.snakemake.log
```

View job-specific logs:
```bash
# Standard output
cat results/logs/{rule}_{wildcards}.out

# Standard error
cat results/logs/{rule}_{wildcards}.err

# Find logs with errors
ls -lt results/logs/*.err | head -10
grep -l 'Error\|error\|Traceback' results/logs/*.err
```

**Finding log files from Snakemake output:**

When Snakemake shows a failed job like:
```
rule geneRegion2Genewise:
    wildcards: seqid=group17400
```

The log file is at: `results/logs/geneRegion2Genewise_seqid=group17400.err`

Pattern: `results/logs/{rule}_{wildcards}.err`

### Expected Runtime (Toy Data)

| Stage | Approximate Time |
|-------|------------------|
| RepeatMasking | 15-30 min |
| RNA-seq alignment | 30-60 min |
| Transcript assembly | 20-40 min |
| Homology search | 30-60 min |
| Augustus training | 1-2 hours |
| Gene model combination | 30-60 min |
| EvidenceModeler | 30-60 min |
| **Total** | **4-8 hours** |

> **Note:** Actual runtime depends on cluster load, resource availability, and data size. The toy dataset (~18.6 Mb) runs significantly faster than real genome annotations.

## 6. Post-Processing with TidyGFF

After the pipeline completes, use TidyGFF to clean and standardize the final GFF3 output:

```bash
singularity exec sylvan.sif TidyGFF -i results/{prefix}/{prefix}.final.gff3 -o {prefix}.tidied.gff3
```

## Troubleshooting

### Common Issues

**1. Job fails with "out of memory"**
- Increase memory allocation in config file for the specific rule
- Recommend **4GB per thread** (e.g., 48 threads = 192g)
- Ensure `ncpus` and `threads` values match for each rule
- Check `results/logs/{rule}_{wildcards}.err` for details

**2. Singularity bind path errors**
- Ensure all input paths are accessible from within the container
- Use absolute paths in configuration

**3. Missing input files**
- Verify Git LFS files are downloaded: `git lfs pull`
- Check file paths in configuration match actual locations

**4. Augustus training fails**
- Ensure sufficient training genes (minimum 500 recommended)
- Check `use_augustus` parameter to skip training with pre-trained species

### Getting Help

- Report issues: https://github.com/plantgenomicslab/Sylvan/issues
- Check existing documentation in the `docs/` directory
