#!/usr/bin/env Rscript

#---- LIBRARIES ----#
library(tidyverse)

#---- ARGUMENTS ----#
args <- commandArgs(trailingOnly = T)
clusters_file <- args[1]
blastn_file <- args[2]
lengths_file <- args[3]

#---- LOAD & JOIN CLUSTER & LENGTH DATA ---- #
clusters <- read.csv(clusters_file, header = F, col.names = c("seq","taxa","segment","cluster")) %>%
  group_by(taxa, segment, cluster) %>%
  count() %>%
  ungroup() %>%
  mutate(seq = paste(taxa,segment,cluster,sep="-")) %>%
  select(seq, taxa, segment, cluster, n)

lengths <- read.csv(lengths_file, header = F, col.names = c("seq","length"))

cluster_lengths <- full_join(clusters, lengths, by = "seq")

#---- LOAD BLASTN RESULTS & ADD MISSING ----#
blastn <- read_tsv(blastn_file, col_names = c("taxa","segment","qaccver", "saccver","qlen","slen","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"))
# get list of pairwise comparisons, add back to blastn, filter by taxa & segment, then calculate ID
q_filt <- cluster_lengths %>%
  mutate(qaccver = seq,
         qtaxa = taxa,
         qseg = segment,
         qlength = length) %>%
  select(qaccver, qtaxa, qseg, qlength)
s_filt <- cluster_lengths %>%
  mutate(saccver = seq,
         staxa=taxa,
         sseg=segment,
         slength = length) %>%
  select(saccver,staxa, sseg, slength)
blastn <- cluster_lengths %>%
  select(seq) %>%
  mutate(qaccver=seq,
         tmp='this works') %>%
  pivot_wider(names_from = "seq", values_from = "tmp" ) %>%
  pivot_longer(names_to = "saccver", values_to = "tmp", 2:ncol(.)) %>%
  select(qaccver, saccver) %>%
  full_join(blastn, by = c("qaccver", "saccver")) %>%
  full_join(q_filt, by = "qaccver") %>%
  full_join(s_filt, by = "saccver") %>%
  filter(qtaxa == staxa & qseg == sseg) %>%
  mutate(taxa = qtaxa,
         segment = qseg,
         qlen = qlength,
         slen = slength) %>%
  drop_na(qaccver, saccver) %>%
  select(-qtaxa, -qseg, -qlength, -staxa, -sseg) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  group_by(qaccver, saccver, taxa, segment, qlen) %>%
  summarise(expected_length = max(slen), 
            aligned = sum(length), 
            unaligned = abs(expected_length - aligned), 
            mismatch = sum(mismatch), 
            gapopen = sum(gapopen), 
            pident = 100*(expected_length - (unaligned+mismatch+gapopen)) / expected_length 
          ) %>%
  ungroup()
write.csv(x=blastn, file = "full-blast.csv", quote = F, row.names = F)
#---- PLOT MATRIX ----#
plot_matrix <- function(ts){
  # subset datafralsme & create plot
  df <- blastn %>%
    mutate(taxa_seg = paste(taxa,segment, sep = "-")) %>%
    filter(taxa_seg == ts)
  p <-ggplot(df, aes(x=qaccver, y=saccver, fill = pident))+
      geom_tile()+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90))
  
  n <- df$qaccver %>% unique() %>% length()
  if(n > 10){
    dims <- n*0.5
  }else{
    dims <- 5
  }
  ggsave(p, filename = paste0(ts,".jpg"), dpi = 300, height = dims, width = dims)
}

ts_list <- blastn %>%
  mutate(taxa_seg = paste(taxa,segment, sep = "-")) %>%
  .$taxa_seg %>%
  unique()
dev_null <- lapply(ts_list, FUN = plot_matrix)

#---- JOIN FINAL DATASETS & SUMMARIZE ----#
blastn <- blastn %>%
  filter(qaccver != saccver) %>%
  group_by(qaccver, taxa, segment) %>%
  summarise(min_pident = round(min(pident), digits = 1), max_pident = round(max(pident), digits = 1)) %>%
  rename(seq = qaccver) %>%
  ungroup() %>%
  select(-taxa, -segment)

clusters %>%
  full_join(blastn, by = "seq") %>%
  full_join(lengths, by = "seq") %>%
  select(seq,taxa,segment,cluster,n,length,min_pident,max_pident) %>%
  write.csv(file = "summary.csv", quote = F, row.names = F)