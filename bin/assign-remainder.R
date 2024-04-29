#!/usr/bin/env Rscript

library(tidyverse)

#---- ARGUMENTS ----#
args <- commandArgs(trailingOnly = T)
reps_path <- args[1]
remainder_mash_path <- args[2]
remainder_list_path <- args[3]
taxa_name <- args[4]
segment_name <- args[5]

# set output file name
file.name <- paste(taxa_name,segment_name,"assigned",sep="-")

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
  filter(dist == max(dist)) %>%
  left_join(df.reps, by = "rep") %>%
  select(seq,cluster) %>%
  mutate(taxa = taxa_name,
         segment = segment_name) %>%
  select(seq, taxa, segment, cluster)

#---- DETERMINE SEQS THAT DID NOT ASSIGN ----#
df.not_assigned <- df.all %>%
  filter(!(seq %in% df.assigned$seq))

# if only one sequence was not assigned then add it as its own cluster
if(nrow(df.not_assigned) == 1){
    df.assigned <- df.not_assigned %>%
      mutate(taxa = taxa_name,
             segment = segment_name,
             cluster = min(df.assigned$cluster)+1) %>%
             rbind(df.assigned) %>%
  select(seq, taxa, segment, cluster)
      df.not_assigned <- data.frame()
}

#---- WRITE RESULTS ----#
write.csv(x=df.assigned, file = paste0(file.name,".csv"), row.names = F, quote = F)
write.csv(x=df.not_assigned, file = "not-assigned.csv", row.names = F, quote = F, col.names = F)

  
