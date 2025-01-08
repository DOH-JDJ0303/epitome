#!/usr/bin/env Rscript
version <- "1.0"

# assign-remainder.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

library(tidyverse)

#---- ARGUMENTS ----#
args <- commandArgs(trailingOnly = T)
reps_path <- args[1]
remainder_mash_path <- args[2]
remainder_list_path <- args[3]
taxon_name <- args[4]
segment_name <- args[5]

#---- VERSION ----#
if(reps_path == "version"){
  cat(version, sep = "\n")
  quit(status=0)
}

# set output file name
file.name <- paste(str_replace_all(taxon_name, pattern = ' ', replacement = '_'),segment_name,"assigned",sep="-")

#---- LOAD INPUTS ----#
df.reps <- read.csv(reps_path, header=F) %>%
  rename(rep=1,
         cluster=2)
df.mash <- read.csv(remainder_mash_path, header=F) %>%
  rename(rep=1,
         seq=2,
         dist=3)
df.all <- read.csv(remainder_list_path,header=F) %>%
  rename(seq=1)

#---- ASSIGN SEQS TO CLUSTERS ----#
df.assigned <- df.mash %>%
  group_by(seq) %>%
  filter(dist == min(dist)) %>%
  left_join(df.reps, by = "rep") %>%
  select(seq,cluster) %>%
  mutate(taxon = taxon_name,
         segment = segment_name) %>%
  select(seq, taxon, segment, cluster)

#---- DETERMINE SEQS THAT DID NOT ASSIGN ----#
df.not_assigned <- df.all %>%
  filter(!(seq %in% df.assigned$seq))

# if only one sequence was not assigned then add it as its own cluster
if(nrow(df.not_assigned) == 1){
    df.assigned <- df.not_assigned %>%
      mutate(taxon = taxon_name,
             segment = segment_name,
             cluster = min(df.assigned$cluster)+1) %>%
             rbind(df.assigned) %>%
  select(seq, taxon, segment, cluster)
      df.not_assigned <- data.frame()
}

#---- WRITE RESULTS ----#
write.csv(x=df.assigned, file = paste0(file.name,".csv"), row.names = F, quote = F)
write.csv(x=df.not_assigned, file = "not-assigned.csv", row.names = F, quote = F, col.names = F)

  
