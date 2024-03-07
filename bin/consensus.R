#!/usr/bin/env Rscript

library(tidyverse)
library(Biostrings)

args <- commandArgs(trailingOnly = T)
prefix <- args[1]
aln_path <- args[2]

consensus <- readDNAStringSet(filepath = aln_path) %>%
  data.frame() %>%
  dplyr::rename(seq=1) %>%
  unique() %>%
  group_by(seq) %>%
  separate(col = 1, into = paste0("site_",1:nchar(.)), sep = 1:nchar(.)) %>%
  pivot_longer(names_to = "site", values_to = "base", cols = 1:ncol(.)) %>%
  group_by(site, base) %>%
  count() %>%
  ungroup() %>%
  group_by(site) %>%
  filter(n == max(n)) %>%
  dplyr::slice(1) %>%
  filter(base == "A" | base == "T" | base == "C" | base == "G") %>%
  mutate(site = as.numeric(str_remove_all(site, pattern = 'site_' ))) %>%
  arrange(site) %>%
  drop_na() %>%
  .$base %>%
  paste(collapse="")

filename <- paste0(prefix,".fa")
write(paste0(">",prefix,"\n",consensus), file = filename)
write(nchar(consensus), file = "length.csv")