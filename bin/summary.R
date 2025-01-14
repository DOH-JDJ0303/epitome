#!/usr/bin/env Rscript
version <- "2.0"

# summary.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#---- ARGUMENTS ----#
args <- commandArgs(trailingOnly = T)
input_qc_file <- args[1]
clusters_file <- args[2]
refs_file     <- args[3]
metadata_file <- args[4]

#---- VERSION ----#
if(input_qc_file == "version"){
  cat(version, sep = "\n")
  quit(status=0)
}

#---- LIBRARIES ----#
library(tidyverse)

#----- FUNCTIONS -----#
# function for collapsing rows by column
collapseCols <- function(val){
  res <- val %>%
    str_remove_all(pattern = '[\\[\\]]') %>%
    str_split(pattern = '; ') %>%
    unlist() %>%
    data.frame() %>%
    rename(var = 1) %>%
    filter(!(var %in% c('','NA','none','null','NULL'))) %>%
    drop_na(var) %>%
    unique() %>%
    .$var
  if( length(res) == 0 ){
    res <- ''
  }else if( length(res) == 1 ){
    res <- unlist(res)
  }else if( any(!is.na(suppressWarnings(as.numeric(res)))) ){
    # cat(paste0('\n',as.character(res),' is numeric!\n'))
    res <- paste(min(as.numeric(res)),max(as.numeric(res)),sep=' - ')
  }else{
    # cat(paste0('\n',as.character(res),' is a character!\n'))
    res <- paste0('[',paste(res, collapse = '; '),']')
  }
  return(res)
}
# function for counting number of elements in a bracket list
countElements <- function(data){
    return(str_count(data, pattern = ';') + 1)
}
# function for expanding accessions associated with each QC'd sequence
expandAccessions <- function(row, df){
  df[row,]$accessions %>%
    str_remove_all(pattern = "[\\[\\]]") %>%
    str_split(pattern = '; ') %>%
    data.frame() %>%
    rename(accession = 1) %>%
    mutate(seq = df[row,]$seq) %>%
    return()
}
#---- LOAD DATA ---- #
# required data
df.qc       <- read_csv(input_qc_file)
df.clusters <- read_csv(clusters_file)
df.refs     <- read_csv(refs_file)
# optional data
if(file.exists(metadata_file)){
  df.meta <- read_csv(metadata_file)
}else{
  df.meta <- df.qc %>% 
    select(seq)
}
#----- CLEAN DATA -----#
df.qc <- df.qc %>%
  select(seq, accessions, length) %>%
  rename(input_lengths = length) %>%
  mutate_all(as.character) %>%
  drop_na(seq)
df.clusters <- df.clusters %>%
  select(-taxon, -segment) %>% 
  mutate_all(as.character) %>%
  drop_na(cluster)
df.refs <- df.refs %>%
  rename(n_qc = n) %>%
  mutate_all(as.character) %>%
  drop_na(cluster)
df.meta <- df.meta %>%
  mutate(across(everything(), ~ gsub("[\r\n]", "", .))) %>%
  mutate_all(as.character) %>%
  drop_na(accession) %>%
  filter(accession != 'null') %>%
  select(-taxon, -segment) %>%
  unique()
# expand accessions associated with each QC'd sequence
df.accessionKey <- do.call(rbind, lapply(1:nrow(df.qc), FUN = expandAccessions, df.qc)) %>%
  mutate_all(as.character) %>%
  drop_na(accession) %>%
  unique()

#----- MERGE & CLEAN MORE -----#
df.summary <- df.accessionKey %>%
  full_join(df.meta, by = 'accession') %>%
  drop_na(accession) %>%
  select(-accession) %>%
  full_join(df.qc, by = 'seq') %>%
  drop_na(seq) %>%
  full_join(df.clusters, by = 'seq') %>%
  drop_na(seq) %>%
  full_join(df.refs, by = "cluster") %>%
  mutate(cluster = case_when(is.na(cluster) & ! is.na(seq) ~ 'failed_qc',
                             TRUE ~ cluster),
         ref     = case_when(cluster == 'failed_qc' ~ 'failed_qc',
                             TRUE ~ ref),
         seq = paste0('seq',seq)) %>%
  drop_na(cluster) %>%
  group_by(cluster) %>%
  unique() %>%
  mutate_at(vars(-group_cols()), collapseCols) %>%
  unique() %>%
  ungroup() %>%
  mutate(n_raw = countElements(accessions),
         seq = str_remove_all(seq, pattern = 'seq'))
# fix any missing taxon and segment info
taxonName <- df.summary %>%
  filter(!(segment %in% c('','NA','none','null','NULL'))) %>%
  drop_na(taxon) %>%
  slice(1) %>%
  .$taxon
segmentName <- df.summary %>%
  filter(!(segment %in% c('','NA','none','null','NULL'))) %>%
  drop_na(segment) %>%
  slice(1) %>%
  .$segment
df.summary <- df.summary %>%
  mutate(taxon   = taxonName,
         segment = segmentName) %>%
  drop_na(taxon)
# rename columns
df.summary <- df.summary %>%
  rename(assembly = ref)
# reorder columns
main_cols <- c('taxon','segment','assembly','length','n_qc','n_raw','cluster','condensed')
extra_cols <- df.summary %>%
  select(-any_of(main_cols)) %>%
  colnames()
df.summary <- df.summary %>%
  select(all_of(c(main_cols, extra_cols))) %>%
  select_if(function(col) !all(is.na(col) | col == ""))

#----- WRITE TO FILE -----#
# full summary
write.csv(x = df.summary, file = 'summary.full.csv', quote = T, row.names = F)
# simple summary
main_cols <- c('taxon','segment','assembly','length','next_closest_ani')
extra_cols <- df.meta %>%
  select(-any_of(main_cols)) %>%
  colnames()
df.simple <- df.summary %>%
  filter(!(assembly %in% c('condensed','failed_qc'))) %>%
  select(any_of(c(main_cols, extra_cols))) %>%
  unique() %>%
  unique
write.csv(x = df.simple, file = 'summary.simple.csv', quote = T, row.names = F)
# taxon summary
n_raw <- df.summary %>%
  .$n_raw %>%
  as.numeric() %>%
  sum()
n_pass <- df.summary %>%
  filter(!(assembly %in% c('failed_qc'))) %>%
  .$n_raw %>%
  as.numeric() %>%
  sum()
max_ani <- df.summary %>%
  filter(!(assembly %in% c('condensed','failed_qc'))) %>%
  .$next_closest_ani %>%
  as.numeric() %>%
  max()
df.taxon <- data.frame(taxon   = taxonName,
                       segment = segmentName,
                       n_raw   = n_raw,
                       n_pass  = n_pass,
                       max_ani = max_ani) %>%
            unique()
write.csv(x = df.taxon, file = 'summary.taxon.csv', quote = T, row.names = F)