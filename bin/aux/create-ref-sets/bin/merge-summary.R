#!/usr/bin/env Rscript

# merge-summary.R
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
summarizeColumns <- function(val){
  res <- val %>%
    str_remove_all(pattern = '[\\[\\]]') %>%
    str_split(pattern = '; ') %>%
    unlist() %>%
    data.frame() %>%
    rename(var = 1) %>%
    filter(!(var %in% c('','NA','none','null','NULL'))) %>%
    drop_na(var) %>%
    unique() %>%
    .$var
  if( length(res) == 0 ){
    res <- ''
  }else if( length(res) == 1 ){
    res <- unlist(res)
  }else if( any( str_detect(res, pattern = '[0-9]+ - [0-9]+') ) ){
    # cat(paste0('\n',as.character(res),' is numeric!\n'))
    res <- res %>%
      str_split(pattern = ' - ') %>%
      unlist() %>%
      data.frame() %>%
      rename(var = 1) %>%
      filter(!(var %in% c('','NA','none','null','NULL'))) %>%
      drop_na(var) %>%
      unique() %>%
      .$var 
    res <- paste(min(as.numeric(res)),max(as.numeric(res)),sep=' - ')
  }else if( any(!is.na(suppressWarnings(as.numeric(res)))) ){
    # cat(paste0('\n',as.character(res),' is numeric!\n'))
    res <- paste(min(as.numeric(res)),max(as.numeric(res)),sep=' - ')
  }else{
    # cat(paste0('\n',as.character(res),' is a character!\n'))
    res <- paste0('[',paste(res, collapse = '; '),']')
  }
  return(res)
}
accessionLists <- function(t, s, accessions){
  t <- unique(t)
  s <- unique(s)
  df <- accessions %>%
    str_remove_all(pattern = '[\\[\\]]') %>%
    str_split(pattern = '; ') %>%
    unlist() %>%
    data.frame() %>%
    rename(accession = 1) %>%
    mutate(taxon = t,
           segment = s) %>%
    unique()
    
    write.csv(x = df, file = paste0('accessions/',t,'-',s,'.accessions.txt'), quote = T, row.names = F)
}
columnPrefix <- function(columns, selected_columns, prefix) {
  existing_columns <- intersect(columns, selected_columns)
  columns[columns %in% existing_columns] <- paste(prefix, columns[columns %in% existing_columns], sep = "")
  return(columns)
}
totSpecies <- function(val){
  species <- val %>%
    str_remove_all(pattern = '[\\[\\]]') %>%
    str_split(pattern = "; ") %>%
    unlist() %>%
    unique()

  lspecies <- species %>%
    paste(collapse = ', ')
  nspecies <- length(species)

  return(list(lspecies, nspecies))
}

#----- MAIN -----#
cols.main <- c('taxon','assembly','segment', 'n_raw','n_qc','accessions','collectionDate', 'species', 'subtype', 'genotype', 'serotype', 'lineage', 'clade', 'biotype', 'subgroup', 'group', 'type', 'geographicRegion', 'organismName_host')
df.main <- do.call(bind_rows, lapply(list.files(path = "./", recursive = T, pattern = ".csv", full.names = T), FUN = mergeMetadata)) %>%
  filter(assembly != 'failed_qc') %>%
  filter(is.na(condensed)) %>% 
  select(any_of(cols.main)) %>%
  select_if(function(col) !all(is.na(col) | col == ""))

#----- REFERENCE SHEET-----#
cols.refsheet <- c('taxon','assembly','segment', 'species', 'collectionDate', 'geographicRegion', 'organismName_host', 'subtype', 'genotype', 'serotype', 'lineage', 'clade', 'biotype', 'subgroup', 'group', 'type')
df.refsheet <- df.main %>%
  select(any_of(cols.refsheet)) %>%
  drop_na(assembly) %>%
  group_by(assembly) %>%
  mutate_at(vars(-group_cols()), summarizeColumns) %>%
  unique() %>%
  ungroup() %>%
  mutate(assembly = paste0('references/',assembly,'.fa.gz'))
write.csv(x = df.refsheet, file = "refsheet.csv", quote = T, row.names = F)

#----- TOTALS -----#
md.totals <- df.main %>%
  summarize(tot_raw = sum(as.numeric(n_raw), na.rm = T), tot_qc = sum(as.numeric(n_qc), na.rm = T), tot_refs = n(), tot_species_names = totSpecies(species)[1], tot_species = totSpecies(species)[2]) %>%
  kable()
write(x = md.totals, file = 'totals.md')
df.totals.taxon <- df.main %>%
  group_by(taxon) %>%
  summarize("No. Input Sequences" = sum(as.numeric(n_raw), na.rm = T), "No. Sequences After QC" = sum(as.numeric(n_qc), na.rm = T), "No. References" = n(), "No. Species" = totSpecies(species)[2]) %>%
  ungroup()

#----- MARKDOWN SUMMARY -----#
cols.md <- c('taxon','segment','species','collectionDate','subtype', 'genotype', 'serotype', 'lineage', 'clade', 'biotype', 'subgroup', 'group', 'type')
md.summary <- df.main %>%
  select(any_of(cols.md)) %>%
  group_by(taxon,segment) %>%
  mutate(references = n()) %>%
  ungroup() %>%
  group_by(taxon) %>%
  mutate_at(vars(-group_cols()), summarizeColumns) %>%
  unique() %>%
  ungroup() %>%
  left_join(df.totals.taxon, by = "taxon") %>%
  rename(Taxon = taxon, 
         Segments = segment)
write(x = kable(md.summary), file = 'summary.md')
md.summary.simple <- md.summary %>%
  select(Taxon, Segments, `No. References`, `No. Species`,`No. Input Sequences`)
write(x = kable(md.summary.simple), file = 'summary.simple.md')

#----- ACCESSION LISTS -----#
if( any('accessions' %in% colnames(df.main)) ){
  unlink('accessions', recursive = T)
  dir.create('accessions',recursive = T)
  dev_null <- df.main %>%
    select(taxon, segment, accessions) %>%
    drop_na() %>%
    unique() %>%
    group_by(taxon, segment) %>%
    mutate(dev_null = accessionLists(taxon, segment, accessions))
}
