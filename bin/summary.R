#!/usr/bin/env Rscript
version <- "1.0"

# summary.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#---- ARGUMENTS ----#
args <- commandArgs(trailingOnly = T)
summary_file <- args[1]
fastani_ava_file <- args[2]
fastani_seeds_file <- args[3]
seeds_file <- args[4]
timestamp <- args[5]

#---- VERSION ----#
if(summary_file == "version"){
  cat(version, sep = "\n")
  quit(status=0)
}

#---- LIBRARIES ----#
library(tidyverse)

#---- FUNCTIONS ----#
basename_fa <- function(path){
    result <- basename(path) %>%
      str_remove_all(pattern = ".fa.gz$")
    return(result)
    
}

#---- LOAD SUMMARIES ---- #
clusters <- read.csv(summary_file, header = F, col.names = c("seq","taxa","segment","cluster","n","n2","length")) %>%
  mutate(condensed = case_when(n2 > 1 ~ TRUE,
                               n2 == 1 ~ FALSE)
                               ) %>%
  select(-n2)

#---- LOAD FASTANI AVA RESULTS & ADD MISSING ----#
fastani_ava <- read_tsv(fastani_ava_file, col_names = c("query","ref","ani","mapped","total")) %>%
  select(query, ref, ani) %>%
  mutate(query = basename_fa(query),
         ref = basename_fa(ref))
# get list of pairwise comparisons, add back to fastani_ava, filter by taxa & segment, then calculate ID
q_filt <- clusters %>%
  select(seq, taxa, segment, length) %>%
  rename(query = seq)
r_filt <- clusters %>%
  mutate(ref = seq,
         staxa=taxa,
         sseg=segment) %>%
  select(ref, staxa, sseg)
fastani_ava <- clusters %>%
  select(seq) %>%
  mutate(query=seq,
         tmp='this works') %>%
  pivot_wider(names_from = "seq", values_from = "tmp" ) %>%
  pivot_longer(names_to = "ref", values_to = "tmp", 2:ncol(.)) %>%
  select(query, ref) %>%
  full_join(fastani_ava, by = c("query", "ref")) %>%
  full_join(q_filt, by = "query") %>%
  full_join(r_filt, by = "ref") %>%
  filter(taxa == staxa & segment == sseg) %>%
  drop_na(query, ref) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  select(-staxa, -sseg)
write.csv(x=fastani_ava, file = "full-ani.csv", quote = F, row.names = F)
#---- PLOT MATRIX ----#
plot_matrix <- function(ts){
  # subset datafralsme & create plot
  df <- fastani_ava %>%
    mutate(taxa_seg = paste(taxa,segment, sep = "-")) %>%
    filter(taxa_seg == ts)
  p <-ggplot(df, aes(x=query, y=ref, fill = ani))+
      geom_tile()+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90))
  
  n <- df$query %>% unique() %>% length()
  if(n > 10){
    dims <- n*0.5
  }else{
    dims <- 5
  }
  if(n<100){
    ggsave(p, filename = paste0(ts,".jpg"), dpi = 300, height = dims, width = dims, limitsize = FALSE)
  }else(cat("Matrix image not saved. Too many references!\n"))
  
}

ts_list <- fastani_ava %>%
  mutate(taxa_seg = paste(taxa,segment, sep = "-")) %>%
  .$taxa_seg %>%
  unique()
dev_null <- lapply(ts_list, FUN = plot_matrix)

#---- SUMMARIZE BLAST AVA FURTHER ----#
fastani_ava <- fastani_ava %>%
  filter(query != ref) %>%
  group_by(query, taxa, segment) %>%
  summarise(min_ani = round(min(ani), digits = 1), max_ani = round(max(ani), digits = 1)) %>%
  mutate(min_ani = case_when(min_ani < 80 ~ '< 80',
                             TRUE ~ as.character(min_ani)),
         max_ani = case_when(max_ani < 80 ~ '< 80',
                             TRUE ~ as.character(max_ani))) %>%
  rename(seq = query) %>%
  ungroup() %>%
  select(-taxa, -segment)

#---- JOIN FINAL DATASETS & SUMMARIZE ----#
## WITHOUT SEEDS
clusters %>%
  full_join(fastani_ava, by = "seq") %>%
  select(seq,taxa,segment,cluster,n,condensed,length,min_ani,max_ani) %>%
  write.csv(file = paste0(timestamp,"-summary.csv"), quote = F, row.names = F)

## WITH SEEDS
if(file.exists(fastani_seeds_file) & file.exists(seeds_file)){
  seeds <- read.csv(seeds_file) %>% 
    rename(seed = 1,
           ref = 2) %>%
    mutate(ref = basename_fa(ref))
  read_tsv(fastani_seeds_file, col_names = c("query","ref","ani","mapped","total")) %>%
    select(query, ref, ani) %>%
    mutate(query = basename_fa(query),
           ref = basename_fa(ref)) %>%
    group_by(query) %>%
    filter(ani == max(ani)) %>%
    ungroup() %>%
    filter(ani >= 95) %>%
    mutate(ani = round(ani, digits = 1)) %>%
    rename(seq = query,
           seed_ani = ani) %>%
    full_join(clusters, by = "seq") %>%
    left_join(seeds, by = "ref") %>%
    select(seq,taxa,segment,cluster,n,length,min_ani,max_ani, seed, seed_ani) %>%
    write.csv(file = paste0(timestamp,"-summary.csv"), quote = F, row.names = F)
}