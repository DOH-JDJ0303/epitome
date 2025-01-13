#!/usr/bin/env Rscript

version <- "1.0"

# input-qc.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#---- ARGUMENTS ----#
args <- commandArgs(trailingOnly=T)
taxon       <- args[1]
segment     <- args[2]
seqs_path   <- args[3]
amb_thresh  <- as.numeric(args[4])
len_thresh  <- as.numeric(args[5])
max_cluster <- as.numeric(args[6])

#---- VERSION ----#
if(taxon == "version"){
  cat(version, sep = "\n")
  quit(status=0)
}

#---- LIBRARY ----#
library(tidyverse)
library(Biostrings)

#----- MAIN -----#
# File basename
file.basename <- paste(str_replace_all(taxon, pattern = ' ', replacement = '_'), segment,sep='-')

# Load sequences
seqs <- readDNAStringSet(seqs_path)
# Combine into a data frame 
df.seqs <- data.frame(accession = str_remove(names(seqs), pattern = "\\s.*"), length = width(seqs), seqString = as.character(seqs, use.names = F), stringsAsFactors = FALSE) %>%
  mutate(seqString = toupper(seqString))

# Consolidate sequences
df.seqs <- df.seqs %>%
  group_by(seqString, length) %>%
  reframe(accessions = paste0('[',paste(accession, collapse = '; '),']')) %>%
  ungroup()
# Gather metrics
legalBases <- '[-ATCGRYSWKMBDHVN]'
median.length <- median(df.seqs$length)
filterTest <- function(test){
    if(test){
        return('fail')
    }else {
       return('pass')
    }
}
df.seqs <- df.seqs %>%
  mutate(seq = row_number()) %>%
  group_by(seqString) %>%
  mutate(illegalBases       = str_remove_all(seqString, pattern = legalBases ),
         ambRatio           = str_count(seqString, 'N') / length,
         filter_illegalBases = nchar(illegalBases) > 0, 
         filter_ambRatio    = ambRatio > amb_thresh, 
         filter_length      = sum(length >= median.length*(1+len_thresh) || length <= median.length*(1-len_thresh)) > 0,
         filters            = paste0('[illegalBases: ',filterTest(filter_illegalBases),', ambRatio: ',filterTest(filter_ambRatio),', lenght: ',filterTest(filter_length),']')
         ) %>%
  ungroup()


# Save passing sequences to file
## function for saving fasta file
saveFasta <- function(df,name){
    fasta <- df %>%
      mutate(record = paste0('>',seq,'\n',seqString)) %>%
      .$record %>%
      paste(collapse = '\n')
    write(fasta, name)
}
## filter to only passing sequences & save
df.seqs.pass <- df.seqs %>%
  filter(!(filter_illegalBases) & !(filter_ambRatio) & !(filter_length))
saveFasta(df.seqs.pass, paste0(file.basename,'.all.fa'))
## check if the number of passing sequences exceeds the maximum number of sequences allowed in a cluster
if( nrow(df.seqs.pass) > max_cluster ){
    ## randomize dataframe
    set.seed(11)
    df.seqs.pass <- df.seqs.pass[sample(nrow(df.seqs.pass)),]
    ## get top set & save
    df.top <- df.seqs.pass[1:max_cluster,]
    saveFasta(df.top, paste0(file.basename,'.top.fa'))
    ## get remainder & save
    df.remainder <- df.seqs.pass[(max_cluster+1):nrow(df.seqs.pass),]
    saveFasta(df.remainder, paste0(file.basename,'.remainder.fa'))
}

# Save metadata file
metadata <- df.seqs #%>%
  # select(-filter_illegalBases, -filter_ambRatio, -filter_length)
write.csv(x = metadata, file = paste0(file.basename,'.qc.csv'), quote = T, row.names = T)
