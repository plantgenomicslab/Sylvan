# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Sylvan is a two-phase genome annotation pipeline for plant and other eukaryotic genomes. It integrates multiple evidence sources (RNA-seq, protein homology, ab initio prediction, liftover) into consensus gene models, then applies semi-supervised random forest filtering to produce a high-quality annotation.

## Environment Setup

- **Interactive shell**: `micromamba activate sylvan`
- **Non-interactive (Claude Code Bash / CI)**: `micromamba run -n sylvan <command>` or use the env Python absolute path `/Users/wyim/micromamba/envs/sylvan/bin/python`

## Common Commands

### Run the annotation pipeline (Phase 1)
```bash
./bin/annotate.sh
```

### Run the filter pipeline (Phase 2)
```bash
./bin/filter.sh
```

### Run with toy data (A. thaliana chr4)
```bash
./bin/annotate_toydata.sh
./bin/filter_toydata.sh
```

### Run unit tests
```bash
python bin/test_feature_importance.py
```

### Run pylint (matches CI configuration)
```bash
pylint $(find bin -name "*.py" -type f) \
  --disable=C0114,C0115,C0116 \
  --disable=R0913,R0914,R0915 \
  --disable=W0612,W0613 \
  --max-line-length=120 \
  --exit-zero
```

### Force rerun of a pipeline
```bash
./bin/annotate.sh --forceall
```

## Architecture

### Two-Phase Pipeline

**Phase 1 — Annotate** (`bin/Snakefile_annotate`, ~2900 lines, 40+ rules):
Generates gene models from multiple evidence types, combined via EVM (Evidence Modeler):
- Repeat masking (RepeatMasker/Modeler)
- RNA-seq alignment and assembly (STAR, HiSat2, StringTie, PsiCLASS, PASA)
- Protein alignment (Miniprot, GeneWise, GMAP)
- Ab initio prediction (Helixer, Augustus)
- Liftover from neighbor species (LiftOff)
- GETA pipeline integration (Augustus training, gene model combination)
- EVM consensus merging
- Output: `results/complete_draft.gff3`

**Phase 2 — Filter** (`bin/Snakefile_filter`, ~460 lines, 20 rules):
Semi-supervised random forest filtering of draft gene models:
- PfamScan (protein domains), RSEM (expression), BLAST (homology), lncDC (lncRNA), BUSCO (conservation)
- Iterative RF training with data-driven heuristic for initial labels
- Output: `results/FILTER/filter.gff3`

### Key Python Modules (all in `bin/`)

| File | Purpose |
|------|---------|
| `Filter.py` | Semi-supervised RF training loop (`semiSupRandomForest`), evidence integration |
| `filter_feature_importance.py` | Leave-one-feature-out ablation analysis |
| `CombineDuck.py` | DuckDB-based gene annotation database for model selection |
| `combine_genemodel.py` | Merges Augustus, TransFrag, GeneWise models with intron support validation |
| `TidyGFF.py` | GFF3 validation, chromosome renaming, output formatting |
| `score_filter.py` | Alternative scoring approach (`bin/Snakefile_filter_score`) |
| `Pick_Primaries.py` | Primary isoform selection |
| `repeat_gene_removal.py` | Repeat-based gene filtering |
| `gff_to_evm.py` | GFF3 to EVM format conversion |

### Configuration

- `config/config_annotate.yml` — All pipeline inputs, parameters, **and** SLURM cluster settings (genome, RNA-seq, proteins, neighbor species, tool settings, per-rule resource allocation). Serves as both pipeline config and `--cluster-config` for Snakemake (like `config_filter.yml`). Replace `placeholder` values before running.
- `config/config_filter.yml` — Filter phase configuration (also serves as `--cluster-config`)
- `config/evm_weights.txt` — Evidence weights for EVM consensus

### Environment Variables

| Variable | Purpose |
|----------|---------|
| `SYLVAN_CONFIG` | Override annotate config path (default: `config/config_annotate.yml`) |
| `SYLVAN_FILTER_CONFIG` | Override filter config path |
| `SYLVAN_RESULTS_DIR` | Override results directory |
| `TMPDIR` | Set automatically to `$(pwd)/results/TMP` by entry scripts |

### Execution Model

- Snakemake 7 orchestrates both phases
- All bioinformatics tools run inside a Singularity container (`singularity/Singularity.def`)
- SLURM cluster execution with per-rule resource allocation
- Entry scripts (`bin/annotate.sh`, `bin/filter.sh`) set up TMPDIR, configure cluster submission, and pass through additional snakemake arguments via `"$@"`

## Languages and Dependencies

- **Python 3.10+** (primary): pandas, scikit-learn, duckdb, pysam, biopython, PyYAML
- **Perl**: `bin/fillingEndsOfGeneModels.pl` (gene model end-filling)
- **R**: `bin/filter_distributions.R` (distribution visualization)
- **Bash**: Entry-point scripts and SLURM job wrappers

## CI

GitHub Actions runs pylint on Python 3.10 and 3.11 for all `bin/*.py` files on pushes/PRs to `main` and `develop`.
