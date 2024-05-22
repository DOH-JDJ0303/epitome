#!/usr/bin/env Rscript
version <- "1.0"

# condense.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#---- ARGUMENTS ----#
args <- commandArgs(trailingOnly = T)
dist_path <- args[1]
lengths_path <- args[2]
clusters_path <- args[3]
taxa_name <- args[4]
segment_name <- args[5]
threshold <- args[6]

#---- VERSION ----#
if(dist_path == "version"){
  cat(version, sep = "\n")
  quit(status=0)
}

#---- LIBRARIES ----#
library(tidyverse)
library(ggtree)
library(ape)

file.base <- paste(taxa_name,segment_name,sep="-")

#---- LOAD CLUSTER SET & GET COUNT ----#
clusters.df <- read_csv(clusters_path) %>%
  mutate(seq = paste(taxa,segment,cluster,sep = "-")) %>%
  group_by(seq,taxa,segment,cluster) %>%
  count()

#---- LOAD SEQ LENGTHS ----#
len.df <- read_csv(lengths_path, col_names = c("seq","length"))

#---- CONDENSE CLOSELY RELATED REFERENCES ----#
# This is only performed when two or more references were generated
if(nrow(clusters.df) > 1){
  #---- LOAD PAIRWISE DISTANCES ----#
  dist.df <- read_tsv(dist_path, col_names = c("ID1","ID2","DIST", "PVAL", "HASH")) %>%
    select(ID1, ID2, DIST)
  dist.mat <- dist.df %>%
    pivot_wider(names_from="ID2", values_from="DIST") %>%
    column_to_rownames(var="ID1") %>%
    as.matrix() %>%
    as.dist()

  #---- HIERARCHICAL CLUSTER ----#
  # create tree
  tree <- hclust(dist.mat, method = "complete") %>%
    as.phylo()
  # test that tree is ultrametric and rooted (required for cutree)
  if(! is.ultrametric(tree)){ 
    cat("Error: Tree is not ultrameric!")
    q()
  }
  if(! is.rooted(tree)){ 
    cat("Error: Tree is not rooted!")
    q()
  }

  #---- CUT DENDROGRAM AT DISTANCE THRESHOLD ----#
  clusters.refs.df <- cutree(as.hclust(tree), h = as.numeric(threshold)) %>%
    data.frame() %>%
    rownames_to_column(var = "seq") %>%
    rename(cluster2 = 2) %>%
    left_join(clusters.df, by = "seq") %>%
    left_join(len.df, by = "seq") %>%
    group_by(cluster2) %>%
    mutate(n2 = n()) %>%
    filter(!(n < 10 & n2 > 1)) %>%
    filter(length == max(length)) %>%
    filter(n2 == max(n2)) %>%
    slice(1) %>%
    ungroup() %>%
    select(seq,taxa,segment,cluster,n,n2,length)
}else{
  clusters.refs.df <- clusters.df %>%
    left_join(len.df, by = "seq") %>%
    mutate(n2 = 1) %>%
    select(seq,taxa,segment,cluster,n,n2,length)
}
#----- SAVE OUTPUT -----#
write.csv(x= clusters.refs.df, file = paste0(file.base,".condensed.csv"), quote = F, row.names = F)



