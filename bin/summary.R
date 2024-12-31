#!/usr/bin/env Rscript
version <- "1.0"

# summary.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#---- ARGUMENTS ----#
args <- commandArgs(trailingOnly = T)
cluster_refs_file <- args[1]
cluster_all_file <- args[2]
raw_seqs_file <- args[3]
clean_seqs_file <- args[4]
metadata_file <- args[5]
fastani_ava_file <- args[6]
prefix <- args[7]

#---- VERSION ----#
if(cluster_refs_file == "version"){
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

#---- LOAD SUMMARY ---- #
df.clusters_refs <- read_csv(cluster_refs_file)

#---- CLASSIFY INPUT SEQS ----#
# function for collapsing rows by group per column
collapse_column <- function(col, df){
    res <- df[,c("cluster",col)] %>%
      rename(original = 2) %>%
      group_by(cluster) %>%
      mutate(collapsed = case_when( is.numeric(original) ~ paste0(min(original)," - ",max(original)),
                                    TRUE ~ paste(unique(original), collapse = "; "))) %>%
      select(cluster, collapsed) %>%
      unique() %>%
      ungroup()
    colnames(res) <- c("cluster",col)
    return(res)
}
# load raw sequences
df.raw_seqs <- read_tsv(raw_seqs_file, col_names = c("id","bases")) %>%
  mutate_all(as.character)
# load cleaned sequences
df.clean_seqs <- read_tsv(clean_seqs_file, col_names = c("seq","bases")) %>%
  mutate_all(as.character)
# merge cluster info with clean sequences and then raw sequences
df.meta <- read_csv(cluster_all_file) %>%
  mutate_all(as.character) %>%
  select(seq, cluster) %>%
  full_join(df.clean_seqs, by = "seq") %>%
  full_join(df.raw_seqs, by = "bases") %>%
  drop_na(seq) %>%
  mutate(cluster = case_when(is.na(cluster) ~ 'filtered sequences',
                             TRUE ~ cluster))
if(file.exists(metadata_file)){
  # load additional metdata & join
  df.meta <- read_tsv(metadata_file) %>%
    rename(id = 1) %>%
    right_join(df.meta, by = "id")
}
meta.cols <- df.meta %>%
  select(-seq, -cluster, -bases) %>%
  colnames()
df.meta <- lapply(meta.cols, FUN = collapse_column, df.meta) %>%
  reduce(left_join, by = "cluster")

#---- LOAD FASTANI AVA RESULTS & ADD MISSING ----#
fastani_ava <- read_tsv(fastani_ava_file, col_names = c("query","ref","ani","mapped","total")) %>%
  select(query, ref, ani) %>%
  mutate(query = basename_fa(query),
         ref = basename_fa(ref))
# get list of pairwise comparisons, add back to fastani_ava, filter by taxa & segment, then calculate ID
q_filt <- df.clusters_refs %>%
  select(seq, taxa, segment, length) %>%
  rename(query = seq)
r_filt <- df.clusters_refs %>%
  mutate(ref = seq,
         staxa=taxa,
         sseg=segment) %>%
  select(ref, staxa, sseg)
fastani_ava <- df.clusters_refs %>%
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

#---- SUMMARIZE FASTANI AVA FURTHER ----#
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
# Join data sources
df.summary <- df.clusters_refs %>%
  mutate_all(as.character) %>%
  full_join(fastani_ava, by = "seq") %>%
  full_join(df.meta, by = "cluster")
# update condensed
taxa_val <- df.summary %>%
  drop_na(taxa) %>%
  slice(1) %>%
  .$taxa
segment_val <- df.summary %>%
  drop_na(segment) %>%
  slice(1) %>%
  .$segment
df.summary <- df.summary %>%
  group_by(cluster) %>%
  mutate(taxa = taxa_val,
         segment = segment_val,
         seq = case_when( is.na(seq) & ! is.na(cluster) ~ 'Condensed' ,
                          TRUE ~ seq ),
         n = case_when( seq == 'Condensed' ~ as.character(length(unlist(str_split(id, pattern = '; ')))),
                        TRUE ~ n )) %>%
  ungroup()
write.csv(x = df.summary, file = paste0(prefix,"-summary.csv"), quote = T, row.names = F)