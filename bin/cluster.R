#!/usr/bin/env Rscript

library(tidyverse)
library(ggtree)

args <- commandArgs(trailingOnly = T)
dist_path <- args[1]
threshold <- args[2]


dist.mat <- read_tsv(dist_path, col_names = c("ID1","ID2","DIST","PVAL","HASHES")) %>%
  select(ID1, ID2, DIST) %>%
  pivot_wider(names_from="ID2", values_from="DIST") %>%
  column_to_rownames(var="ID1") %>%
  as.matrix() %>%
  as.dist()

tree <- hclust(dist.mat, method = "complete")

clusters <- cutree(tree, h = as.numeric(threshold)) %>%
  data.frame() %>%
  rownames_to_column(var = "seq") %>%
  rename(cluster = 2)

write.csv(clusters, file = 'clusters.csv', quote = F, row.names = F)

p <- ggtree(tree)%<+%clusters+
  geom_tippoint(aes(color = as.character(cluster)))+
  labs(color = "Cluster")
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
ggsave(plot = p, filename = "clusters.jpg", dpi = 300, width = wdth, height = hght, limitsize=F)