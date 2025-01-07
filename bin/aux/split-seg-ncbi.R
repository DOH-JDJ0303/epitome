library(tidyverse)
library(Biostrings)

args <- commandArgs(trailingOnly=T)
meta_path <- args[1]
seqs_path <- args[2]

meta <- read_csv(meta_path)
seqs <- readDNAStringSet(seqs_path)

filename.base <- basename(seqs_path) %>%
  str_remove_all(pattern = ".fa")

split_fasta <- function(seg){
    meta.seg <- meta %>%
      filter(Segment == seg) %>%
      drop_na(Segment)
    seq_list <- meta.seg$Accession
    seqs.seg <- seqs[seq_list]

    filename.seg <- paste0(filename.base,"_seg-",seg,".fa")
    writeXStringSet(seqs.seg, filename.seg)
}

dev_null <- lapply(unique(meta$Segment), FUN = split_fasta)