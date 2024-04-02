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
threshold <- args[4]

#---- VERSION ----#
if(dist_path == "version"){
  cat(version, sep = "\n")
  quit(status=0)
}

#---- LIBRARY ----#
library(tidyverse)
library(ggtree)
library(ape)

# set output file name
file.name <- paste(taxa_name,segment_name,sep="-")

#---- LOAD PAIRWISE DISTANCES ----#
dist.mat <- read_tsv(dist_path, col_names = c("ID1","ID2","DIST","PVAL","HASHES")) %>%
  select(ID1, ID2, DIST) %>%
  pivot_wider(names_from="ID2", values_from="DIST") %>%
  column_to_rownames(var="ID1") %>%
  as.matrix() %>%
  as.dist()

#---- HIERARCHICAL CLUSTER ----#
# create tree
tree <- hclust(dist.mat, method = "average") %>%
  as.phylo()
# test that tree is ultrametric (required for cutree)
if(! is.ultrametric(tree)){ 
  cat("Error: Tree is not ultrameric!")
  q()
}
# convert back to hclust for cutree
tree <- as.hclust(tree)
#---- CUT DENDROGRAM AT DISTANCE THRESHOLD ----#
clusters <- cutree(tree, h = as.numeric(threshold)) %>%
  data.frame() %>%
  rownames_to_column(var = "seq") %>%
  rename(cluster = 2) %>%
  mutate(taxa = taxa_name,
         segment = segment_name) %>%
  select(seq, taxa, segment, cluster)

# adjust height to percentage for more intuitive interpretation
tree$height <- 100*tree$height
x_scale <- seq(from = -max(round(tree$height+0.5, digits = 0)), to = 0, by = 10)

# plot tree
p <- ggtree(tree)%<+%clusters+
  geom_tippoint(aes(color = as.character(cluster)))+
  geom_vline(xintercept = -100*as.numeric(threshold), linetype = "dashed", color = "#E35335")+
  theme_tree2()+
  labs(color = "Cluster", x = "Approx. Nucleotide Difference (%)")+
  scale_x_continuous(breaks=x_scale, labels=abs(x_scale))
#---- SAVE OUTPUTS ----#
# cluster info
write.csv(clusters, file = paste0(file.name,'.csv'), quote = F, row.names = F)
n_iso <- p$data %>%
  drop_na() %>%
  nrow()

# set image dimensions
wdth <- n_iso/100
if(wdth < 200){
  wdth <- 10
}
hght=n_iso/100
if(hght < 200){
  hght <- 10
}
# plot
ggsave(plot = p, filename = paste0(file.name,'.jpg'), dpi = 300, width = wdth, height = hght, limitsize=F)