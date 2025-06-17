##############################################################
## Pauline Trinh April 25, 2024 
## Purpose: Check the gisaid files for influenza B, filter out coinfections
## pulled all fluB sequences from Jan 1, 2014 to Apr 23, 2024 
## Refer to the following s3 bucket for files/folders used: 
## s3://wamep-nf-bucket/viral_references/
##############################################################

##############################################################
## Sequences were downloaded manually 
## based on submission dates. Submission dates for each .fasta time period
## have been written into the names of the fasta files. 
## original, complete (all 8 segments sequenced) sequences were selected 
## Download Sequences (DNA) as FASTA using the following header scheme: 
## Isolate ID | DNA Accession no. | Isolate name | Segment | Segment number
##############################################################
##############################################################
library(ape)
library(tidyverse)

fluB_metadata <- readxl::read_excel("gisaid_fluB_2014_2024/raw_gisaid_pull/2014_Jan1_2024_Apr23_gisaid_fluB_metadata.xls")
dim(fluB_metadata)
#17481 x 64

fluB_fasta_header_data <- read.FASTA("gisaid_fluB_2014_2024/raw_gisaid_pull/2014_Jan1_2024_Apr23_gisaid_fluB_sequence.fasta") %>% 
  labels() %>% 
  tibble() %>%
  rename_at(1, ~'header') %>% 
  separate(col=header, 
           into = c("Isolate_Id","dna_accession","name","segment", "segment_number"), 
           sep="\\|", 
           remove = FALSE) 
dim(fluB_fasta_header_data) #140436      6

# let's check how many segments were sequenced per isolate_id 
check_fluB_segments <- fluB_fasta_header_data %>% 
  group_by(Isolate_Id) %>% 
  summarise(count_segments = n()) %>% 
  arrange(desc(count_segments))
# IN total there are 17481 isolates of fluB
# 80 isolates seem to be coinfections so we'll filter those out of the fasta file and metadata 

### Filter out coinfection isolates from metadata 
no_coinfections_fluB <- fluB_fasta_header_data %>% 
  group_by(Isolate_Id) %>% 
  summarise(count_segments = n()) %>% 
  filter(count_segments == 8)
# 17401 isolates remaining 

fluB_seqs_withmetadata <- no_coinfections_fluB %>% 
  left_join(fluB_metadata)

write.table(fluB_seqs_withmetadata, "gisaid_fluB_2014_2024/2014_Jan1_2024_Apr23_gisaid_fluB_metadata.txt")
write.csv(fluB_seqs_withmetadata, "gisaid_fluB_2014_2024/2014_Jan1_2024_Apr23_gisaid_fluB_metadata.csv", row.names = F)
saveRDS(fluB_seqs_withmetadata, "gisaid_fluB_2014_2024/2014_Jan1_2024_Apr23_gisaid_fluB_metadata.rds")

# generate a .txt file that we'll use seqtk subseq to subsample the isolates without a co-infection 
fluB_seqtk_seqs <- fluB_fasta_header_data %>% 
  full_join(check_fluB_segments) %>% 
  filter(count_segments == 8) %>% 
  select(header)
dim(fluB_seqtk_seqs) #139208 sequences that we'll need to use corresponding to 17401 isolates 
write.table(fluB_seqtk_seqs, "gisaid_fluB_2014_2024/fluB_seqtk_seqs.txt", sep= "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)



####################################
## Double check cleaned data files 
####################################

# check_fasta_file_headers <- read.FASTA("gisaid_fluB_2014_2024/2014_Jan1_2024_Apr23_gisaid_fluB_sequences_cleaned.fasta") %>% 
#   labels() %>% 
#   tibble() %>%
#   rename_at(1, ~'header') %>% 
#   separate(col=header, 
#            into = c("Isolate_Id","dna_accession","name","segment", "segment_number"), 
#            sep="\\|", 
#            remove = FALSE) 
# dim(check_fasta_file_headers)# should be 139208 sequences. check!
# 
# check_metadata_file_txt <- read.table("gisaid_fluB_2014_2024/2014_Jan1_2024_Apr23_gisaid_fluB_metadata.txt")
# check_metadata_file_csv <- read.csv("gisaid_fluB_2014_2024/2014_Jan1_2024_Apr23_gisaid_fluB_metadata.csv")
# check_metadata_file_rds <- readRDS("gisaid_fluB_2014_2024/2014_Jan1_2024_Apr23_gisaid_fluB_metadata.rds")
