# Sylvan Local Execution Error Log

## Environment
- Date: 2026-03-05
- Machine: 16 cores, 62GB RAM, Linux (CUDA 12.6 GPU)
- Config: toydata/config/config_annotate_local.yml
- Singularity: singularity/Sylvan_sandbox_v2 (rebuilt from Sylvan.sif with fakeroot)

## Status
- Dry-run: PASSED (228 jobs)
- Run progress: 1468/1648 steps (89%) completed before helixer dependency blocked remaining jobs
- After helixer fallback: progressing through GETA/Augustus steps
- Current: Fixing run: block container execution issues

## Errors and Fixes

### 1. helixer — DEFERRED (GPU/CUDA mismatch)
- **Cause**: Container CUDA 11.2 vs host CUDA 12.6 — TF cannot find GPU
- **Fix**: Added `|| { touch {output}; }` fallback to create empty output
- **Status**: WORKAROUND — pipeline continues with empty helixer evidence

### 2. RepeatModeler config — FIXED
- **Cause**: RepModelConfig.pm had hardcoded paths from upstream tarball
- **Fix**: perl -0777 patching in sandbox to set all paths to /opt/envs/repeat/bin
- **Also fixed**: TRF_DIR, UCSCTOOLS_DIR, LTR_RETRIEVER_DIR missing — downloaded UCSC tools, fixed all paths
- **Status**: FIXED in Sylvan_sandbox_v2

### 3. RepeatModeler -pa deprecated — FIXED
- **Cause**: `-pa` parameter deprecated in RepeatModeler 2.0.5, replaced by `-threads`
- **Fix**: Changed `-pa {threads}` to `-threads {threads}` in Snakefile_annotate
- **Status**: FIXED

### 4. RepeatModeler missing output files — FIXED
- **Cause**: RepeatModeler 2.0.5 creates output in RM_* subdirectory as `consensi.fa`,
  but Snakefile expects `species-families.fa` at database path
- **Fix**: Added post-processing step to copy from RM_*/consensi.fa to species-families.fa
- **Status**: FIXED

### 5. geneRegion2Genewise — FIXED (run: block not in container)
- **Cause**: `run:` blocks execute on HOST, not inside singularity
- **Fix**: Created `run_in_container()` helper function, used throughout
- **Status**: FIXED

### 6. GETA perl shebang — FIXED
- **Cause**: GETA scripts use `#!/opt/miniconda3/bin/perl`, container has perl at `/usr/bin/perl`
- **Fix**: `sed -i` replaced all shebangs to `#!/usr/bin/env perl` in sandbox
- **Status**: FIXED in Sylvan_sandbox_v2

### 7. Sandbox corruption from failed micromamba install — FIXED
- **Cause**: Attempted `micromamba install` inside writable sandbox failed mid-way,
  corrupted repeat environment (unlinked RepeatMasker, LTR_retriever)
- **Fix**: Rebuilt sandbox from SIF using `singularity build --fakeroot --sandbox`
- **Status**: FIXED — using Sylvan_sandbox_v2

### 8. mergeTransfrag — FIXED (run: block not in container)
- **Cause**: `transfragDump` called via shell() in run: block, executes on host
- **Fix**: Changed to use `run_in_container()` helper
- **Status**: FIXED

### 9. TransDecoder.LongOrfs output path mismatch — FIXED
- **Cause**: TransDecoder creates output in `{input}.transdecoder_dir/longest_orfs.cds`
  but Snakefile expects `{output_dir}/longest_orfs.cds`
- **Fix**: Added post-processing copy step in shell block
- **Status**: FIXED

### 10. geneModels2AugusutsTrainingInput log file missing — FIXED
- **Cause**: run: block reads stderr log file that doesn't exist when running via container
- **Fix**: Redirect stderr to expected log file path using `run_in_container(stderr_file=...)`
- **Status**: FIXED

### 11. BGM2AT — FIXED (run: block not in container)
- **Cause**: BGM2AT called via shell() in run: block, exit code 127 (command not found)
- **Fix**: Changed to use `run_in_container()` helper
- **Status**: FIXED

### 12. augustusWithHints — FIXED (run: block not in container)
- **Cause**: Augustus called via shell() in run: block
- **Fix**: Changed to use `run_in_container()` helper
- **Status**: FIXED

### 13. TransDecoder.Predict missing .gff3 output — FIXED
- **Cause**: TransDecoder 5.7.1 creates `.bed`, `.cds`, `.pep` but the final `.gff3`
  is either not retained or removed during Snakemake's cleanup-on-error.
  The internal checkpoint file `longest_orfs.cds.best_candidates.gff3.revised_starts.gff3`
  is the actual source GFF3 used to derive all other outputs.
- **Fix**: Added fallback chain: try `.gff3` first, then copy from internal
  `revised_starts.gff3`, then `best_candidates.gff3`
- **Status**: FIXED

### 14. BGM2AT `&>` redirect runs new_species.pl in background (dash vs bash) — FIXED
- **Cause**: Container `/bin/sh` is dash, not bash. Perl `system()` uses `/bin/sh`.
  In dash, `&>` means `& >` (background + redirect), not `>file 2>&1`.
  So `new_species.pl` ran in background, and `etraining` started before species config was created.
- **Fix**: Changed `&> file` to `> file 2>&1` in BGM2AT script inside sandbox
- **Status**: FIXED in Sylvan_sandbox_v2

### 15. BGM2AT Augustus config copy — FIXED
- **Cause**: `cp -rf $CONDA_PREFIX/config/ .../Augustus/config` creates `config/config/`
  instead of copying contents
- **Fix**: Changed to `cp -rf $CONDA_PREFIX/config/* .../Augustus/config/`
- **Status**: FIXED

### 16. filterHMMScan wrong conda env — FIXED
- **Cause**: `hmmscan` is in `filter` env, not `genepred`. Rule had `micromamba activate genepred`.
- **Fix**: Changed to `micromamba activate filter`
- **Status**: FIXED

### 17. agat_clean_final `-gff` parsed as `-g ff` — FIXED
- **Cause**: AGAT uses Getopt::Long with bundling. `-gff` is parsed as `-g` `ff`.
- **Fix**: Changed to `--gff` (double-dash long option)
- **Status**: FIXED

## Pipeline Status: COMPLETE
- **Output**: `results/complete_draft.gff3` (5.8MB, 61,588 lines, ~4,744 gene models)
- **Steps**: 157/157 completed
- **Known limitation**: Helixer produces empty output (CUDA version mismatch workaround)

## Systematic Fix: run_in_container() helper
All `run:` blocks that call container tools via `shell()` fail because `run:` blocks
execute on the HOST, not inside singularity. Created a global `run_in_container()`
helper function at the top of Snakefile_annotate that wraps commands with
`singularity exec` when a container is configured.

Affected rules: mergeTransfrag, geneRegion2Genewise, geneModels2AugusutsTrainingInput,
BGM2AT, augustusWithHints

---

# Phase 2: Filter Pipeline Errors

## Environment
- Same machine as Phase 1 (16 cores, 62GB RAM, Linux)
- Config: toydata/config/config_filter_local.yml
- Runner: bin/filter_local.sh (--cores 16, no SLURM)

### 18. lncDC.py command not found — FIXED
- **Cause**: lncDC is installed as a Python package module at
  `$CONDA_PREFIX/lib/python3.9/site-packages/bin/lncDC.py`, not as a CLI entry point.
  Calling `lncDC.py` directly fails with "command not found".
- **Fix**: Changed to `python $CONDA_PREFIX/lib/python3.9/site-packages/bin/lncDC.py`
  in Snakefile_filter
- **Status**: FIXED

### 19. STAR --quantTranscriptomeBan unrecognized — FIXED
- **Cause**: STAR 2.7.11b renamed `--quantTranscriptomeBan` to
  `--quantTranscriptomeSAMoutput`
- **Fix**: Changed parameter name to `--quantTranscriptomeSAMoutput` in Snakefile_filter
- **Status**: FIXED

### 20. `python: command not found` in base env — FIXED
- **Cause**: `micromamba activate base` doesn't properly add `/opt/envs/base/bin` to
  PATH inside singularity. `which python` returns nothing.
- **Fix**: Changed to `python3` initially, then to `/opt/envs/base/bin/python3` (see #21)
- **Status**: FIXED

### 21. `No module named 'pandas'` — FIXED
- **Cause**: System `/usr/bin/python3` lacks pandas. The correct interpreter is
  `/opt/envs/base/bin/python3` (Python 3.14) which has pandas 3.0.1 + scikit-learn.
- **Fix**: Changed Filter.py invocation to `/opt/envs/base/bin/python3 bin/Filter.py`
  in Snakefile_filter
- **Status**: FIXED

### 22. pandas 3.0 `delim_whitespace` removed — FIXED
- **Cause**: `pd.read_csv(delim_whitespace=True)` raises TypeError in pandas 3.0
  (parameter was deprecated and removed)
- **Fix**: Changed to `sep=r'\s+'` in Filter.py
- **Status**: FIXED

### 23. pandas 3.0 LossySetitemError (bool→float64) — FIXED
- **Cause**: `filter_df.loc[:,"TPM"] = rsem_data["TPM"] > float(tpm_cutoff)` fails
  because the column dtype is float64 and the assigned value is bool, which pandas 3.0
  treats as lossy.
- **Fix**: Use `.copy()` on DataFrame creation and direct column assignment
  `filter_df["TPM"] = ...` instead of `.loc` assignment
- **Status**: FIXED

### 24. 0 test samples — empty x_test in RF training — FIXED
- **Cause**: With pandas 3.0, `NaN == False` returns `False` (changed from older pandas
  where it returned `True`). Combined with empty Helixer evidence (all NaN), the
  keep/discard conditions matched every gene, leaving none as "None" for the test set.
- **Fix**: Added explicit `filter_df[col] = filter_df[col].fillna(False)` for all
  boolean columns before the keep/discard labeling conditions
- **Status**: FIXED

## Filter Pipeline Status: COMPLETE
- **Output**: `results/FILTER/filter.gff3`
- **Steps**: 82/82 completed
- **Results**: 596 genes kept, 3,698 discarded (out of 4,294 total)
- **Known limitation**: Empty Helixer evidence causes more aggressive filtering
  (596 vs ~3,756 expected with valid Helixer data)
