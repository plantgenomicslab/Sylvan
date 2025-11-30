library(ggplot2)
library(readr)
library(dplyr)

setwd("/Users/jslomas/Box/Yim_Lab/GeneModelFilter_old")

rex <- read_tsv("BLASTP.OUT.Rex", col_names = FALSE, col_select=c(1,3,13)) %>% 
  setNames(c("gene_id", "rex_pident", "rex_qcov"))
blast <- read_tsv("BLASTP.OUT.TMP", col_names = FALSE, col_select=c(1,3,13)) %>% 
  setNames(c("gene_id", "blast_pident","blast_qcov"))
rsem <- read_tsv("RSEM.genes.results", col_names = TRUE, col_select=c("gene_id", "TPM"))
augustus <- read_tsv("augustus_coverage.bed", col_names = FALSE, col_select=c(5,9))%>% 
  setNames(c("gene_id", "augustus"))
helixer <- read_tsv("helixer_coverage.bed", col_names = FALSE, col_select=c(5,9))%>% 
  setNames(c("gene_id", "helixer"))
rsem_cov <- read_tsv("cov.bed", col_names = FALSE, col_select=c(1,7))%>% 
  setNames(c("gene_id", "rsem_cov"))
repeat_cov <- read_tsv("repeat_coverage.bed", col_names = FALSE, col_select=c(5,9))%>% 
  setNames(c("gene_id", "repeat_cov"))

ggplot(rex, aes(x=rex_pident)) + geom_histogram()
ggplot(rex, aes(x=rex_qcov)) + geom_histogram()
ggplot(blast, aes(x=blast_pident)) + geom_histogram()
ggplot(blast, aes(x=blast_qcov)) + geom_histogram()
ggplot(rsem, aes(x=TPM)) + geom_histogram() + scale_y_log10() + scale_x_log10()
ggplot(augustus, aes(x=augustus)) + geom_histogram()
ggplot(helixer, aes(x=helixer)) + geom_histogram()
ggplot(rsem_cov, aes(x=rsem_cov)) + geom_histogram()
ggplot(repeat_cov, aes(x=repeat_cov)) + geom_histogram()
