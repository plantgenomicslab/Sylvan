#!/bin/bash

# This script is intened to be run after both the annotation phase 
# and filter phases have been completed. It will clean up intermediate 
# files generated during the annotation phase. Primary tool outputs
# used for EVM are retained. The logs and configs are also preserved 
# for later reference. 


echo "Removing untracked annotation phase workflow files."
snakemake --list-untracked -s bin/Snakefile_annotate 2>&1 | \
  grep -v LncDC_plant_db | \
  grep -v bin | \
  grep -v config | \
  grep -v singularity | \
  grep -v README.md | \
  grep -v LICENSE | \
  grep -v FILTER | \
  grep -v complete_draft | \
  grep -v logs | \
  xargs -r rm -f

find . -type d -empty -delete
rm -rf .RepeatMaskerCache/ pasaPost.ok hmm_files_bak pasa.sqlite.alt_splice_label_combinations.dat

echo "Removing trimmed RNAseq files"
shopt -s nullglob

files=(GETA/fastp/single/*.fq.gz GETA/fastp/paired/*.fq.gz)

if [ ${#files[@]} -gt 0 ]; then
  rm -f "${files[@]}"
fi

shopt -u nullglob
