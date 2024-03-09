#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = T)
clusters_file <- args[1]
blastn_file <- args[2]

clusters <- read.csv(clusters_file, header = F, col.names = c("seq","taxa","segment","cluster")) %>%
  group_by(taxa, segment, cluster) %>%
  count() %>%
  ungroup() %>%
  mutate(seq = paste(taxa,segment,cluster,sep="-")) %>%
  select(seq, cluster, n)
blastn <- read_tsv(blastn_file, col_names = c("taxa","segment","qaccver", "saccver","qlen","slen","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")) %>%
  filter(qaccver != saccver) %>%
  group_by(taxa, segment, qaccver) %>%
  summarize(length = max(qlen), min_aln = min(length), max_aln = max(length), min_pident = min(pident), max_pident = max(pident)) %>%
  ungroup() %>%
  rename(seq = qaccver)

clusters %>%
  full_join(blastn, by = "seq") %>%
  select(seq,taxa,segment,cluster,n,length,min_aln,max_aln,min_pident,max_pident) %>%
  write.csv(file = "summary.csv", quote = F, row.names = F)