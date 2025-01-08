#!/usr/bin/env Rscript
version <- "1.1"

# condense.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#---- ARGUMENTS ----#
args <- commandArgs(trailingOnly = T)
dist_path <- args[1]
lengths_path <- args[2]
clusters_path <- args[3]
taxon_name <- args[4]
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

file.base <- paste(str_replace_all(taxon_name, pattern = ' ', replacement = '_'),segment_name,sep="-")

#---- LOAD CLUSTER SET & GET COUNT ----#
clusters.df <- read_csv(clusters_path) %>%
  filter(taxon == taxon_name & segment == segment_name) %>%
  mutate(ref = paste(str_replace_all(taxon,pattern = ' ', replacement = '_'),segment,cluster,sep = "-")) %>%
  group_by(ref,taxon,segment,cluster) %>%
  count()

#---- LOAD SEQ LENGTHS ----#
len.df <- read_csv(lengths_path, col_names = c("ref","length"))

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
    rownames_to_column(var = "ref") %>%
    rename(cluster2 = 2) %>%
    left_join(clusters.df, by = "ref") %>%
    left_join(len.df, by = "ref") %>%
    group_by(cluster2) %>%
    mutate(condensed = case_when(n() > 1 ~ paste0('[',paste(cluster, collapse = ","),']'),
                                 TRUE ~ '[]')) %>%
    ungroup()
  #---- SELECT BEST REFERENCE FOR EACH CLUSTER ----#
  select.refs <- clusters.refs.df %>%
    filter(!(n < 10 & n() > 1)) %>%
    filter(length == max(length)) %>%
    filter(n() == max(n())) %>%
    slice(1) %>%
    ungroup() %>%
    select(ref,taxon,segment,cluster,n,condensed,length)
  #----- NEXT CLOSEST REFERENCE -----#
  # select references - considers only other selected references
  next_closest <- dist.df %>%
    filter(ID1 != ID2) %>%
    filter(ID1 %in% select.refs$ref) %>%
    filter(ID2 %in% select.refs$ref) %>%
    drop_na() %>%
    group_by(ID1) %>%
    filter(DIST == min(DIST,na.rm = T)) %>%
    reframe(next_closest_ani = round(100*(1-DIST)), next_closest_ref = paste0('[',paste(ID2, collapse = '; '),']')) %>%
    ungroup() %>%
    rename(ref = ID1)
  # condensed references - consider select references
  if(nrow(select.refs) != nrow(clusters.refs.df)){
    next_closest <- dist.df %>%
      filter(ID1 != ID2) %>%
      filter(!(ID1 %in% select.refs$ref)) %>%
      filter(ID2 %in% select.refs$ref) %>%
      drop_na() %>%
      group_by(ID1) %>%
      filter(DIST == min(DIST,na.rm = T)) %>%
      reframe(next_closest_ani = round(100*(1-DIST)), next_closest_ref = paste0('[',paste(ID2, collapse = ', '),']')) %>%
      ungroup() %>%
      rename(ref = ID1) %>%
      rbind(next_closest)
  }
  # merge
  clusters.refs.df <- clusters.refs.df %>%
    left_join(next_closest, by = "ref")
  #---- ASSIGN OTHER REFERENCES AS CONDENSED ----#
  clusters.refs.df <- clusters.refs.df %>%
    mutate(ref = case_when(!(ref %in% select.refs$ref) ~ 'condensed',
                           TRUE ~ ref )) %>%
    select(ref,taxon,segment,cluster,n,condensed,length,next_closest_ani,next_closest_ref) %>%
    unique()

}else{
  clusters.refs.df <- clusters.df %>%
    left_join(len.df, by = "ref") %>%
    mutate(condensed = '[]') %>%
    select(ref,taxon,segment,cluster,n,condensed,length)
}

#----- SAVE OUTPUT -----#
write.csv(x= clusters.refs.df, file = paste0(file.base,".condensed.csv"), row.names = F)



