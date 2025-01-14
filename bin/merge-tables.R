#!/usr/bin/env Rscript

# merge-tables.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

library(tidyverse)

mergeMetadata <- function(file){
    df <- read_csv(file) %>%
      mutate(across(everything(), ~ gsub("[\r\n]", "", .))) %>%
      mutate_all(as.character)
    return(df)
}

result <- do.call(bind_rows, lapply(list.files(path = "./", recursive = T, pattern = ".csv", full.names = T), FUN = mergeMetadata))

write.csv(x = result, file = "merged.csv", quote = T, row.names = F)