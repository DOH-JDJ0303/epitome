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

#---- LOAD DATA ---- #
# required data
df.qc       <- read_csv(input_qc_file)
df.clusters <- read_csv(clusters_file)
df.refs     <- read_csv(refs_file)
# optional data
if(file.exists(metadata_file)){
  # expand accessions associated with each QC'd sequence
  expandAccessions <- function(row, df){
    df[row,]$accessions %>%
      str_remove_all(pattern = "[\\[\\]]") %>%
      str_split(pattern = ',') %>%
      data.frame() %>%
      rename(accession = 1) %>%
      mutate(seq = df[row,]$seq) %>%
      return()
  }
  df.accessionKey <- do.call(rbind, lapply(1:nrow(df.qc), FUN = expandAccessions, df.qc))

  # prepare metadata
  df.meta <- read_csv(metadata_file) %>%
    merge(df.accessionKey, by = 'accession') %>%
    drop_na(accession) %>%
    unique() %>%
    ungroup() %>%
    select(-taxon, -segment)
}else{
  df.meta <- data.frame(seq = df.qc$seq)
}

#----- CLEAN DATA -----#
df.qc <- df.qc %>%
  select(seq, accessions, length) %>%
  rename(input_lengths = length) %>% 
  mutate_all(as.character)
df.clusters <- df.clusters %>%
  select(-taxon, -segment) %>% 
  mutate_all(as.character)
df.refs <- df.refs %>%
  rename(n_qc = n) %>%
  mutate_all(as.character)
df.meta <- df.meta %>%
  mutate_all(as.character)


#----- MERGE & CLEAN MORE -----#
collapseCols <- function(data){
    data <- data %>%
      str_split(pattern = ", ") %>%
      unlist() %>%
      str_remove_all(pattern = '[\\[\\]]') %>%
      unique()
    data <- data[grep(pattern = '(null|NA)', data, invert = T)] %>%
      paste(collapse = ', ')
    return(paste0('[', data, ']'))
}
getRange <- function(data){
  data <- str_split(data, pattern = ', ') %>%
    unlist() %>%
    str_remove_all(pattern = '[\\[\\]]') %>%
    as.numeric()
  if( length(data) > 1 ){
    return(paste(min(data),max(data),sep=' - '))
  }else{
    return(as.character(data))
  }
}
cleanBrackets <- function(data){
  return(str_remove_all(data, pattern = '[\\[\\]]'))
}
countElements <- function(data){
  data %>%
    str_remove_all(pattern = '[\\[\\]]') %>%
    str_split(pattern = ', ') %>%
    unlist() %>%
    length() %>%
    return()
}

df.summary <- df.qc %>%
  merge(df.meta, by = "seq") %>%
  drop_na(seq) %>%
  merge(df.clusters, by = 'seq') %>%
  drop_na(seq) %>%
  merge(df.refs, by = "cluster") %>%
  group_by(cluster) %>%
  unique() %>%
  mutate_all(collapseCols) %>%
  unique() %>%
  mutate(input_lengths = getRange(input_lengths),
         taxon         = cleanBrackets(taxon),
         segment       = cleanBrackets(segment),
         length        = cleanBrackets(length),
         n_qc          = cleanBrackets(length),
         cluster       = cleanBrackets(cluster),
         ref           = cleanBrackets(ref),
         n_raw         = countElements(accessions) ) %>%
  drop_na(taxon)
# reorder columns
main_cols <- c('taxon','segment','ref','length','n_qc', 'n_raw', 'cluster', 'condensed')
extra_cols <- df.summary %>%
  select(-main_cols) %>%
  colnames()
df.summary <- df.summary %>%
  select(main_cols, extra_cols)

#----- WRITE TO FILE -----#
write.csv(x = df.summary, file = 'summary.csv', quote = T, row.names = F)

nrow(df.accessionKey)
sum(df.summary$n_raw)




