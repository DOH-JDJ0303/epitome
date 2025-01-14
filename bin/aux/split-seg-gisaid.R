#!/usr/bin/env Rscript

# split-seg-gisaid.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#---- ARGUMENTS ----#
args <- commandArgs(trailingOnly=T)
fasta_file    <- args[1]
metadata_file <- args[2]
taxon         <- args[3]
species       <- args[4]

#---- LIBRARY ----#
library(tidyverse)
library(Biostrings)

#----- MAIN -----#
# fasta
seqs <- readDNAStringSet(fasta_file)
names(seqs) <- data.frame(name = names(seqs)) %>%
  group_by(name) %>%
  mutate(name = as.character(unlist(str_split(name, pattern = '\\|'))[2])) %>%
  .$name
# metadata
meta <- read_csv(metadata_file)
filename.base <- basename(metadata_file) %>%
  str_remove(pattern = '\\.csv')
segmentKey <- data.frame(segment_name = c("PB2","PB1","PA","HA","NP","NA","MP","NS"), segment = 1:8)
select_cols <- c('taxon', 'geographicRegion','collectionDate','organismName_host','subtype','clade','lineage','HA_type','NA_type')
meta <- meta %>%
  group_by(Isolate_Id) %>%
  mutate(taxon             = taxon,
         species           = species,
         geographicRegion  = unlist(str_split(Location, ' / '))[1],
         collectionDate    = substr(Collection_Date, 1, 4),
         organismName_host = Host,
         serotype          = str_remove_all(Subtype, pattern = 'A( / )'),
         lineage           = Lineage,
         clade             = Clade,
         HA_type           = str_extract(subtype, pattern = 'H[1-9]+'),
         NA_type           = str_extract(subtype, pattern = 'N[1-9]+') ) %>%
  ungroup() %>%
  select(all_of(c(paste0(segmentKey$segment_name, ' Segment_Id'),select_cols))) %>%
  pivot_longer(names_to = 'segment_name', values_to = 'accession', 1:8) %>%
  group_by(accession) %>%
  mutate(segment_name = str_remove_all(segment_name, pattern = ' Segment_Id'),
         accession    = unlist(str_split(accession, pattern = '\\|'))[1]) %>%
  ungroup() %>%
  left_join(segmentKey, by = 'segment_name') %>%
  filter(accession %in% names(seqs)) %>%
  drop_na(segment)
# save results
saveBySeg <- function(seg){
    meta.seg <- meta %>%
      filter(segment == seg) %>%
      drop_na(accession)
    seg_name <- unique(meta.seg$segment_name)
    seg_num  <- unique(meta.seg$segment)
    filename.seg <- paste0(filename.base,"_seg-",seg_num,"-",seg_name)
    writeXStringSet(seqs[meta.seg$accession], paste0(filename.seg,'.fa'))
    meta.seg <- meta.seg %>%
      select(any_of(c('accession','segment',select_cols)))
    write.csv(x = meta.seg, file = paste0(filename.seg,'.csv'), quote = T, row.names = F)
}
dev_null <- lapply(unique(meta$segment), FUN = saveBySeg)


