# Sylvan Genome Annotation

Sylvan computes comprehensive genome annotations by combining and filtering the results of EVM/PASA, GETA, Mikado, and Helixer.
Sylvan is designed to run via the SLURM job scheduler on a cloud or other high-performance computing system.

Sylvan runs in two phases, the annotation phase and the filtering phase.
The annotation phase processes all the provided sources of evidence (RNAseq, protein databases, neighbor species, and *ab initio* predictions), computes multiple draft genome annotations, and merges them into a single draft.
The filtering phase uses a semi-supervised random forest algorithm to remove spurious gene models.

The pipeline also includes **TidyGFF**, a robust script for formatting GFF files according to user defined naming conventions.
After running the pipeline, use `bin/TidyGFF.py` to ready your GFF for public distribution.

## Table of Contents
- [Installation and environment](#installation-and-environment)
  - [Uploading to Sylabs Cloud (for developers)](#uploading-to-sylabs-cloud-for-developers)
- [Toy data](#toy-data)
- [Preparing for a run](#preparing-for-a-run)
  - [Input data download](#input-data-download)
  - [Construct a repeat library](#construct-a-repeat-library)
  - [Configure the annotation and filter phases](#configure-the-annotation-and-filter-phases)
- [Running the pipeline](#running-the-pipeline)
  - [Important: Singularity directory binding](#important-singularity-directory-binding)
  - [SLURM configuration](#slurm-configuration)
  - [Dry run and execution](#dry-run-and-execution)
- [Logs and troubleshooting](#logs-and-troubleshooting)
- [Formatting output with TidyGFF](#formatting-output-with-tidygff)

## Installation and environment
---------------
The quickest way to obtain Sylvan is to download the singularity image, clone the repo, and install Snakemake with conda.
```bash
# Download the singularity image from Sylabs Cloud Library
singularity pull --arch amd64 sylvan.sif library://wyim/sylvan/sylvan:latest

# Download the repo
git clone https://github.com/plantgenomicslab/Sylvan.git
cd Sylvan

# Create an environment with Snakemake version 7
conda create -n snake -c conda-forge -c bioconda \
  python=3.11 \
  snakemake=7 -y

conda activate snake
```

If you have root permission and you prefer to build the singularity image from source:
```bash
# Build the singularity image
cd Sylvan/singularity

sudo singularity build \
  --bind ~/Sylvan/singularity/.bashrc:/.bashrc \
  Sylvan.sif Sylvan.def
```

### Uploading to Sylabs Cloud (for developers)

To upload a rebuilt container to [Sylabs Cloud Library](https://cloud.sylabs.io):

```bash
# Login to Sylabs Cloud (create access token at https://cloud.sylabs.io/tokens)
singularity remote login SylabsCloud

# Push the container (use -U for unsigned push)
singularity push -U sylvan.sif library://wyim/sylvan/sylvan:latest
```

> **Note:** The `-U` flag pushes unsigned containers. For signed containers, see the [Sylabs documentation](https://docs.sylabs.io/guides/latest/user-guide/signNverify.html).

## Toy data

A toy dataset is included in the `toydata/` directory for testing and demonstration purposes. The dataset contains:

- **genome.fasta**: *A. thaliana* chromosome 4 split into 3 equal segments (~18.6 Mb total)
- **Neighbor species**: Chr4 data from *A. lyrata* (aly4), *C. rubella* (cru4), *C. hirsuta* (chi4), and *A. thaliana* (ath4)
- **neighbor.cds**: Combined CDS sequences from syntenic regions
- **Pre-built STAR index**: For RNA-seq mapping

> **_NOTE:_** Toy data files are stored with [Git LFS](https://git-lfs.github.com/). Install Git LFS before cloning, or run `git lfs pull` after cloning to download the files.

For detailed instructions on using the toy data, see the [Wiki](Wiki.md).

## Preparing for a run
-----------------

### Input data download

**Required inputs for the annotation phase:**

| Input | Description | Config Field |
|-------|-------------|--------------|
| Genome assembly | FASTA file (`.fa`, `.fasta`, `.fa.gz`, or `.fasta.gz`) | `genome` |
| RNA-seq data | Gzipped FASTQ files in a single folder | `rna_seq` |
| Protein sequences | FASTA file(s) from UniProt, OrthoDB, etc. | `proteins` |
| Neighbor species | GFF3 annotations and genome assemblies | `liftoff.neighbor_gff`, `liftoff.neighbor_fasta` |
| Repeat library | Repeat nucleotide FASTA file (e.g., from EDTA) | `RM_lib` |
| Singularity image | Path to `sylvan.sif` | `singularity` |

**Additional inputs for the filter phase:**

| Input | Description | Config Field (in config_filter.yml) |
|-------|-------------|-------------------------------------|
| RexDB | Plant repeat element database | `RexDB` |
| PFAM database | HMM database for protein domains (included in Singularity) | `HmmDB` |

---

#### 1. Download RNAseq data (fastq - paired/unpaired/both)

All fastq files (all library types) should be placed into a single folder with `.fastq.gz` suffix.

**Supported naming conventions for paired-end reads:**
- `sample_1.fastq.gz` / `sample_2.fastq.gz`
- `sample_R1.fastq.gz` / `sample_R2.fastq.gz`

**Single-end reads:**
- `sample.fastq.gz` (no `_1`, `_2`, `_R1`, or `_R2` suffix)

Symlinks may be used to avoid storing multiple copies of the raw data.

**Example directory structure:**
```
rna_seq/
├── SRR12345_1.fastq.gz
├── SRR12345_2.fastq.gz
├── SRR67890_R1.fastq.gz
├── SRR67890_R2.fastq.gz
└── single_end_sample.fastq.gz
```

#### 2. Download proteins (fasta)

Download protein fasta files from UniProt, OrthoDB, etc. For best results, exclude proteins from closely related species to avoid circular evidence.

#### 3. Download neighbor species annotations (fasta, gff3)

Place all assembly (fasta) files into one folder and all annotation files (GFF3) into another. We recommend organizing these into `neighbors/gff` and `neighbors/genome`. If possible, obtain the CDS sequences as well for use with EDTA, otherwise you will need to generate these yourself.

> **_NOTE:_** Files from the same species must start with a unique three letter code (e.g. `Ath.gff3`, `Ath.fasta`)

**Supported file extensions:**
- **Genome assemblies:** `.fa` or `.fasta` (NOT `.fas`)
- **Annotations:** `.gff` or `.gff3`

**Example directory structure:**
```
neighbors/
├── gff/
│   ├── Ath.gff3
│   ├── Osa.gff3
│   └── Zma.gff3
└── genome/
    ├── Ath.fasta
    ├── Osa.fasta
    └── Zma.fasta
```

#### 4. Download RexDB (for filter phase only)

The filter phase requires a library of plant-specific known repeat elements from RexDB. This is configured in `config_filter.yml`, not `config_annotate.yml`.
```bash
wget https://github.com/repeatexplorer/rexdb/blob/main/Viridiplantae_v4.0.fasta
```

### Construct a repeat library

We use EDTA to construct a repeat sequence annotation for the genome. For EDTA, obtain CDS sequences for the neighbor species you will use in the analysis.

**Download the EDTA container (recommended):**
```bash
export SINGULARITY_CACHEDIR=$PWD
singularity pull EDTA.sif docker://quay.io/biocontainers/edta:2.2.0--hdfd78af_1
```

**Run EDTA:**
```bash
singularity exec --cleanenv --env PYTHONNOUSERSITE=1 EDTA.sif EDTA.pl \
    --genome genome.fasta \
    --cds neighbors.cds \
    --anno 1 \
    --threads 16
```

The output file `genome.fa.mod.EDTA.TElib.fa` should be used as the `RM_lib` parameter in the configuration.

### Configure the annotation and filter phases

Both phases are managed by YAML configuration files. Copy them from `config/` to your working directory:
```bash
cp config/config_annotate.yml config/config_filter.yml .
```

---

### Annotation Phase (`config_annotate.yml`)

#### Required Input Parameters

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `genome` | file | Genome assembly (`.fa`, `.fasta`, `.fa.gz`, or `.fasta.gz`) | `genome.fa.gz` |
| `rna_seq` | directory | Folder containing RNA-seq FASTQ files | `rna_seq/` |
| `proteins` | file | Protein FASTA file(s) | `proteins.fa` |
| `singularity` | file | Singularity image or sandbox | `sylvan.sif` |
| `liftoff.neighbor_gff` | directory | Folder with neighbor GFF3 files | `neighbors/gff/` |
| `liftoff.neighbor_fasta` | directory | Folder with neighbor genome FASTA files | `neighbors/genome/` |
| `helixer_model` | string | Helixer model type | `land_plant`, `vertebrate`, `invertebrate`, or `fungi` |
| `helixer_subseq` | int | Helixer subsequence length | `64152` (plants), `21384` (fungi), `213840` (vertebrates) |
| `RM_species` | string | RepeatMasker species | `Embryophyta` |
| `RM_lib` | file | Repeat nucleotide FASTA file | `repeats.fa` |
| `augustus_species` | string | Augustus species name or custom name for training | `ath`, `tomato`, `my_species` |
| `evm_weights` | file | EVM weights file (provided in `config/`) | `config/evm_weights.txt` |

#### Optional Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `use_augustus` | string | Species identifier for Augustus. When set, `augustus_species` is ignored, HMM files must exist, and training is skipped (saves runtime). | `placeholder` (train new) |
| `augustus_start_from` | string | Species identifier to start Augustus optimization from. Training will use this species' parameter file as starting point (may save runtime). | `placeholder` (start fresh) |
| `num_evm_files` | int | Number of EVM partition splits | `126` |

#### SLURM Cluster Configuration

| Parameter | Description |
|-----------|-------------|
| `account` | SLURM billing/allocation account (not your username) |
| `partition` | SLURM partition |

---

### Filter Phase (`config_filter.yml`)

The filter phase uses outputs from the annotation phase plus additional external inputs.

#### Outputs from Annotation Phase (auto-generated)

| Parameter | Type | Description |
|-----------|------|-------------|
| `genome` | file | `results/GETA/genome.fasta` |
| `anot_gff` | file | `results/complete_draft.gff3` |
| `augustus_gff` | file | `results/GETA/Augustus/augustus.gff3` |
| `helixer_gff` | file | `results/AB_INITIO/Helixer/helixer.gff3` |
| `protein` | file | `results/PROTEIN/merged_rmdup_proteins.fasta` |
| `repeat_gff` | file | Repeat annotation GFF from annotation phase |

#### Additional External Inputs

| Parameter | Type | Description |
|-----------|------|-------------|
| `RexDB` | file | RexDB repeat element database (download required) |
| `HmmDB` | directory | PFAM database directory (included in Singularity at `/usr/local/src`) |
| `busco_lin` | string | BUSCO lineage (e.g., `eudicots_odb10`) |

#### SLURM Cluster Configuration

Same as annotation phase - set `account` and `partition` in the `__default__` section.

---

### Finding Your SLURM Account

To find your SLURM account, partition, and QOS settings:
```bash
sacctmgr show user "$USER" withassoc -nP \
  format=User,Cluster,Account,Partition,QOS,DefaultQOS \
  | column -s '|' -t
```

**Example output:**
```
User   Cluster    Account            Partition          QOS      DefaultQOS
wyim   pronghorn  cpu-s1-pgl-0       cpu-s1-pgl-0       normal
wyim   pronghorn  cpu-s6-test-0      cpu-s6-test-0      normal   normal
wyim   pronghorn  gpu-s2-pgl-0       gpu-s2-core-0      normal
```

| Column | Description |
|--------|-------------|
| User | Your username |
| Cluster | HPC cluster name |
| Account | Billing/allocation account (use this for `account` in config) |
| Partition | Available partition (use this for `partition` in config) |
| QOS | Quality of Service level |
| DefaultQOS | Default QOS for this association |

To see available partitions:
```bash
sinfo
```

> **_IMPORTANT:_** The SLURM `account` is your billing/allocation account, **not** your username!


## Running the pipeline
-----------------

### Important: Singularity directory binding

> **WARNING:** Snakemake with Singularity only mounts the current working directory by default. If your input data (genome, RNA-seq, proteins, etc.) resides outside the Sylvan directory, Singularity will not be able to access it and the pipeline will fail.

**Option 1 (Recommended):** Store or symlink all input data within the Sylvan working directory.
```bash
# Example: create symlinks to your data
ln -s /path/to/genome.fa.gz genome.fa.gz
ln -s /path/to/rna_seq rna_seq
ln -s /path/to/proteins.fa proteins.fa
```

**Option 2:** Set the `SINGULARITY_BIND` environment variable to bind additional directories:
```bash
export SINGULARITY_BIND="/path/to/data1,/path/to/data2"
```

**Option 3:** Modify `bin/annotate.sh` to include singularity bind arguments:
```bash
snakemake -p \
    --use-singularity \
    --singularity-args "--bind /path/to/data" \
    ...
```

### SLURM configuration

> **Note:** The SLURM `account` in the configuration file is your **billing/allocation account**, not your HPC username.

**System requirements:**

Sylvan was developed and tested on a SLURM cluster with **64GB RAM** and **64 CPU** nodes. The default resource allocations in `config_annotate.yml` reflect these specifications. You may need to adjust the memory and CPU settings according to your system's capabilities.

**Important configuration adjustments:**

| Setting | Default | Notes |
|---------|---------|-------|
| Job time limit | 14 days | Adjust `time` in `__default__` section based on your cluster's wall-time limits |
| Memory per rule | Varies | Recommend **4GB per thread** (e.g., 48 threads = 192g) |
| CPUs per rule | Varies | Scale according to your node's CPU count; `ncpus` and `threads` should match |

**Resource allocation guidelines:**
```yaml
# - ncpus: Number of CPUs requested from SLURM (--cpus-per-task)
# - threads: Number of threads used by the rule (should match ncpus)
# - memory: Recommend 4g per thread (e.g., 48 threads = 192g)

example_rule:
  ncpus: 48      # SLURM CPU allocation
  threads: 48   # Actual threads used (should match ncpus)
  memory: 192g  # 4g × 48 threads
```

> **_GPU nodes:_** Some clusters require separate `account` and `partition` settings for GPU jobs (e.g., Helixer). If your cluster has dedicated GPU partitions, uncomment and configure the `account` and `partition` fields under the `helixer:` section in the config file.


To find your SLURM account(s):
```bash
sacctmgr show user $USER withassoc format=account
```

To see available partitions (including GPU partitions):
```bash
sinfo
```

To check your account's resource limits:
```bash
sacctmgr show qos format=name,maxwall,maxsubmit
```

### Dry run and execution

Before running each phase, perform a dry run with Snakemake to ensure that the configuration was done properly and all file paths are correct:
```bash
# Dry run - shows what will be executed without running anything
snakemake -np -s bin/Snakefile_annotate
```

If the dry run completes without errors, run the annotation phase:
```bash
# Run the annotation phase
sbatch -A [account] \
    -p [partition] \
    -c 1 \
    --mem=1g \
    -J annotate \
    -o annotate.out \
    -e annotate.err \
    --wrap="bin/annotate.sh"
```

**Expected outputs from annotation phase:**

All outputs are organized under the `results/` directory:
- `results/complete_draft.gff3` - Combined gene models from EVM/PASA, Helixer, Mikado, and GETA
- `results/AB_INITIO/Helixer/helixer.gff3` - Helixer predictions
- `results/GETA/Augustus/augustus.gff3` - Augustus predictions
- `results/LIFTOVER/LiftOff/liftoff.gff3` - Liftoff annotations from neighbor species
- `results/TRANSCRIPT/PASA/pasa.sqlite.pasa_assemblies.gff3` - PASA transcript assemblies
- `results/EVM.all.gff3` - EVM consensus predictions
- `results/logs/` - SLURM job logs for each rule

**Output directory structure:**
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

Check the filter configuration by performing a dry run and then execute `bin/filter.sh` via SLURM:
```bash
# Dry run for the filter phase
snakemake -np -s bin/Snakefile_filter

# Run the filter phase
sbatch -A [account] \
    -p [partition] \
    -c 1 \
    --mem=4g \
    -J filter \
    -o filter.out \
    -e filter.err \
    --wrap="bin/filter.sh"
```

**Final output:** `results/FILTER/filter.gff3`

### Tuning the filter

After the filtration script has completed, the results may be fine-tuned by adjusting the parameters of the random forest.
To adjust hyperparameters, rerun the python script outside of Snakemake:
```bash
singularity exec sylvan.sif python bin/Filter.py --help
```

## Logs and troubleshooting

Given variations in input data size, the program may terminate prematurely due to memory issues. In the event that the program stops prematurely, first check the SLURM log to identify which Snakemake rule produced the error.
```bash
grep Error annotate.err
grep Error filter.err
```

Log files for each rule are stored in `results/logs/`. Look in `results/logs/[rule]_[wildcards].err` and/or `results/logs/[rule]_[wildcards].out` for potential memory issues. If a rule is crashing due to memory constraints, the memory can be adjusted in the `Cluster Configuration` section of each configuration file.

**Finding log files from Snakemake output:**

When Snakemake reports an error for a specific rule, the log file path follows this pattern:
```
results/logs/{rule}_{wildcards}.err
```

For example, if Snakemake shows:
```
rule geneRegion2Genewise:
    wildcards: seqid=group17400
```

The log file is at: `results/logs/geneRegion2Genewise_seqid=group17400.err`

For rules without wildcards (e.g., `liftoff`), the log is: `results/logs/liftoff_.err`

### Common issues

1. **MissingInputException errors:** If Snakemake reports missing input files (e.g., `results/TRANSCRIPT/spades/okayset/spades.okay.tr`), this typically means an upstream rule failed.
   - Verify RNA-seq files are correctly named (`*_1.fastq.gz`/`*_2.fastq.gz` or `*_R1.fastq.gz`/`*_R2.fastq.gz` for paired-end)
   - Check that all input paths are accessible (see Singularity binding section above)
   - Review logs in `results/logs/` for the specific upstream rule that failed
   - Run `snakemake -np -s bin/Snakefile_annotate` to see which files are expected

2. **Singularity "file not found" errors:** Ensure all input data paths are bound to the container (see "Singularity directory binding" section).

3. **SLURM account errors:** The `account` field requires your billing/allocation account, not your username. Run `sacctmgr show user $USER withassoc format=account` to find your account.

4. **Out of memory errors:** Increase memory allocation in the config file for the specific rule that failed. Check `results/logs/[rule]_[wildcards].err` for memory-related errors.

5. **EDTA TensorFlow errors:** If you see `ImportError: cannot import name 'to_categorical'`, use the official EDTA container:
   ```bash
   singularity pull EDTA.sif docker://quay.io/biocontainers/edta:2.2.0--hdfd78af_1
   ```

### Useful Snakemake commands
```bash
# Dry run up to a specific rule
snakemake -np -s bin/Snakefile_annotate --until <rule>

# Mark files as complete without re-running
snakemake --touch -c 1 -s bin/Snakefile_annotate --until <rule>

# Show summary of all rules and their status
snakemake -n -s bin/Snakefile_annotate --summary

# Unlock the working directory (if Snakemake was interrupted)
bin/annotate.sh --unlock
bin/filter.sh --unlock

# Force re-run of a specific rule
snakemake -s bin/Snakefile_annotate --forcerun <rule>

# Generate a visual DAG of the workflow
snakemake -s bin/Snakefile_annotate --dag | dot -Tpdf > dag.pdf
```

## Formatting output with TidyGFF

After annotation is complete, use TidyGFF to format the output GFF for public distribution:

```bash
singularity exec sylvan.sif python bin/TidyGFF.py \
    [prefix] \
    [input.gff] \
    --out [output_name] \
    --splice-name [t|mRNA|etc.] \
    --justify [digit_number] \
    --sort \
    --chrom-regex [chromosome_regex] \
    --source [source]
```

**Example:**
```bash
singularity exec sylvan.sif python bin/TidyGFF.py \
    Slyc \
    results/FILTER/filter.gff3 \
    --out Slyc_v1.0 \
    --splice-name t \
    --justify 5 \
    --sort \
    --chrom-regex "^Chr" \
    --source Sylvan
```

## No HPC Environment?

If you don't have access to an HPC cluster, you can deploy a SLURM cluster on Google Cloud using the [Google Cloud Cluster Toolkit](https://docs.cloud.google.com/cluster-toolkit/docs/quickstarts/slurm-cluster).

This allows you to:
- Create an on-demand SLURM cluster in the cloud
- Scale resources as needed for your annotation jobs
- Pay only for the compute time you use

## Citation

If you use Sylvan in your research, please cite:

> Sylvan: A comprehensive genome annotation pipeline. *Currently under review.*

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For issues and questions, please open an issue on GitHub: https://github.com/plantgenomicslab/Sylvan/issues
