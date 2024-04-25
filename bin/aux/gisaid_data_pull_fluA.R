##############################################################
## Pauline Trinh April 24, 2024 
## Purpose: Check the gisaid files for Influenza A, filter out coinfections
## pulled all fluA sequences from Jan 1, 2014 to Apr 23, 2024 
## Refer to the following s3 bucket for files/folders used: 
## s3://wamep-nf-bucket/viral_references/
##############################################################

##############################################################
## Sequences in raw_gisaid_pull/ were downloaded manually in chunks < 20k 
## based on submission dates. Submission dates for each .fasta time period
## have been written into the names of the fasta files. 
## original, complete (all 8 segments sequenced) sequences were selected 
## sequence files were then concatenated into one .fasta 
## metadata files are merged together in this R script. 
##############################################################
library(tidyverse)
library(ape)
xls_files <- list.files(path = "gisaid_fluA_2014_2024/raw_gisaid_pull/", pattern = "\\.xls$", full.names = TRUE)

xls_list <- lapply(xls_files, function(file) {
  # Read the CSV file into a data frame
  my_data <- readxl::read_excel(file)
    # Return the modified data frame
  return(my_data)
})

## join the list of dataframes into one and do some variable clean up
xls_combined <- do.call(rbind, xls_list) 
## 78,332 isolates 

#############################################
## Pull in fasta headers to check sequences 
## 2014_Jan1_2024_Apr23_gisaid_epiflu_isolates_concatenated.fasta
## is the concatenated file of fasta files located in raw_gisaid_pull
#############################################

# read in 2014_Jan1_2024_Apr23_gisaid_epiflu_isolates_concatenated.fasta and format the headers 
fasta_header_data <- read.FASTA("gisaid_fluA_2014_2024/2014_Jan1_2024_Apr23_gisaid_epiflu_isolates_concatenated.fasta") %>% 
  labels() %>% 
  tibble() %>%
  rename_at(1, ~'header') %>% 
  separate(col=header, 
           into = c("Isolate_Id","dna_accession","name","segment", "segment_number"), 
           sep="\\|", 
           remove = FALSE) 
dim(fasta_header_data) #634896      6

# let's check how many segments each isolate has 
check_segments <- fasta_header_data %>% 
  group_by(Isolate_Id) %>% 
  summarise(count_segments = n()) %>% 
  arrange(desc(count_segments))
# it looks like there are 378 co-infection isolate samples
# let's filter them out of our metadata and our fasta files 

# if we only went with isolates with no co-infections 
# 78703 isolates 
no_coinfections <- fasta_header_data %>% 
  group_by(Isolate_Id) %>% 
  summarise(count_segments = n()) %>% 
  filter(count_segments == 8)

seqs_withmetadata <- no_coinfections %>% 
  left_join(xls_combined)

write.table(seqs_withmetadata, "gisaid_fluA_2014_2024/2014_Jan1_2024_Apr23_gisaid_fluA_metadata.txt")
write.csv(seqs_withmetadata, "gisaid_fluA_2014_2024/2014_Jan1_2024_Apr23_gisaid_fluA_metadata.csv", row.names = F)
saveRDS(seqs_withmetadata, "gisaid_fluA_2014_2024/2014_Jan1_2024_Apr23_gisaid_fluA_metadata.rds")

# generate a .txt file that we'll use seqtk subseq to subsample the isolates without a co-infection 
seqtk_seqs <- fasta_header_data %>%
  full_join(check_segments) %>% 
  filter(count_segments == 8) %>% 
  select(header)
dim(seqtk_seqs) #629624 sequences that we'll need to use corresponding to 78703 isolates 
write.table(seqtk_seqs, "gisaid_fluA_2014_2024/gisaid_seqtk.txt", sep= "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

