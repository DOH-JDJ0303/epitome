#!/usr/bin/env Rscript

#!/usr/bin/env Rscript
version <- "1.0"

# cluster.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#---- ARGUMENTS ----#
args <- commandArgs(trailingOnly = T)
dist_path <- args[1]
taxa_name <- args[2]
segment_name <- args[3]
iteration <- args[4]
threshold <- args[5]

#---- VERSION ----#
if(dist_path == "version"){
  cat(version, sep = "\n")
  quit(status=0)
}

#---- LIBRARY ----#
library(tidyverse)
library(ggtree)
library(ape)
install.packages("bigmemory")
library(bigmemory)

# set output file name
file.name <- paste(taxa_name,segment_name,iteration,sep="-")

#---- LOAD PAIRWISE DISTANCES ----#
dist.df <- read_tsv(dist_path, col_names = c("ID1","ID2","DIST")) %>%
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
clusters <- cutree(as.hclust(tree), h = as.numeric(threshold)) %>%
  data.frame() %>%
  rownames_to_column(var = "seq") %>%
  rename(cluster = 2) %>%
  mutate(seq = as.numeric(seq),
         taxa = taxa_name,
         segment = segment_name) %>%
  select(seq, taxa, segment, cluster)

# adjust height to percentage for more intuitive interpretation
tree$edge.length <- 100*tree$edge.length

# get initial plot to determine the max tip length from root
p <- ggtree(tree)
x_max <- p$data %>%
    filter(isTip) %>%
    arrange(x) %>%
    slice(1) %>%
    .$x

# plot tree
p <- ggtree(tree)%<+%clusters+
  geom_tippoint(aes(color = as.character(cluster)))+
  geom_vline(xintercept = x_max-(100*as.numeric(threshold)/2), linetype = "dashed", color = "#E35335")+
  theme_tree2()+
  labs(color = "Cluster", x = "Estimated Nucleotide Difference (%)")+
  xlim(0,x_max*1.1)

#---- CHECK FOR DISCREPANCIES ----#
join1 <- clusters %>%
  select(seq,cluster) %>%
  rename(ID1 = seq,
         cluster1 = cluster)
join2 <- clusters %>%
  select(seq,cluster) %>%
  rename(ID2 = seq,
         cluster2 = cluster)
discrep <- dist.df %>%
  mutate(DIST = round(DIST, digits = 2)) %>%
  left_join(join1, by = "ID1") %>%
  left_join(join2, by = "ID2") %>%
  filter(cluster1 == cluster2) %>%
  filter(as.numeric(DIST) > as.numeric(threshold))
if(nrow(discrep)){
  cat("Error: The intra-cluster distance exceeds the supplied threshold", sep ="\n")
  quit(status=1)
}

#---- SAVE OUTPUTS ----#
# cluster info
write.csv(clusters, file = paste0(file.name,'.csv'), quote = F, row.names = F)
n_iso <- p$data %>%
  drop_na() %>%
  nrow()

# set image dimensions
dim <- n_iso/200
if(dim < 10){
  dim <- 10
}
if(dim > 50){
  dim <- 50
}
# plot
ggsave(plot = p, filename = paste0(file.name,'.jpg'), dpi = 300, width = dim, height = dim, limitsize=F)