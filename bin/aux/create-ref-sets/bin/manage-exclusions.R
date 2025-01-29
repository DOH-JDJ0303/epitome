#!/usr/bin/env Rscript

# manage-exclusions.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#----- LIBRARIES -----#
library(tidyverse)
library(knitr)

#----- FUNCTIONS -----#
mergeMetadata <- function(file){
    df <- read_csv(file) %>%
      mutate_all(as.character)
    return(df)
}

#----- MAIN -----#
do.call(bind_rows, lapply(list.files(path = "./", recursive = T, pattern = ".csv", full.names = T), FUN = mergeMetadata)) %>%
  drop_na() %>%
  group_by(taxon, reason) %>%
  summarize(exclusions = n()) %>%
  ungroup() %>%
  select(any_of(c('taxon','exclusions','reason'))) %>%
  kable() %>%
  write(file = 'exclusions.md')

