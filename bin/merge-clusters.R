#!/usr/bin/env Rscript

version <- "1.0"

# merge-clusters.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#---- ARGUMENTS ----#
args <- commandArgs(trailingOnly=T)
top_file       <- args[1]
assigned_file  <- args[2]
looseends_file <- args[3]

#---- VERSION ----#
if(top_file == "version"){
  cat(version, sep = "\n")
  quit(status=0)
}

#---- LIBRARY ----#
library(tidyverse)

#----- LOAD DATA -----#
df.top       <- read_csv(top_file)
if(file.exists(assigned_file)){ 
    df.assigned <- read_csv(assigned_file) 
}
if(file.exists(looseends_file)){ 
    df.looseends <- read_csv(looseends_file) 
}

#----- CLEAN DATA -----#
df.top <- df.top %>%
  mutate_all(as.character) %>%
  mutate(clusterMethod = 'main')
source_list <- list(df.top)
if(file.exists(assigned_file)){ 
  df.assigned <- df.assigned %>%
    mutate_all(as.character) %>%
    mutate(clusterMethod = 'assigned')
  source_list <- list(df.top, df.assigned)
}
if(file.exists(looseends_file)){ 
  df.looseends <- df.looseends %>%
    mutate_all(as.character) %>%
    mutate(clusterMethod = 'loose-end',
           cluster       = as.character(as.numeric(cluster) + max(as.numeric(df.top$cluster)))
           )
    source_list <- list(df.top, df.assigned, df.looseends)
}

#----- COMBINE CLUSTERS -----#
df.all <- do.call(bind_rows, source_list) %>%
  drop_na()

#----- WRITE FILE -----#
write.csv(x = df.all, file = 'clusters.csv', quote = T, row.names = F)