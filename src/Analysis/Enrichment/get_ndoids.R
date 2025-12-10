library(HDO.db)
library(tidyverse)
library(org.Hs.eg.db)

### args
args <- commandArgs(trailingOnly = TRUE)

input_degree  <- args[1]
annotated_degree <- args[2]

### Input
degree_df = read.table(
  input_degree,
  sep="\t",
  header=T
)

### Processing
entrez_id <- mapIds(
  org.Hs.eg.db, keys = degree_df$gene_name,
  column = "ENTREZID", keytype = "SYMBOL")

gene_df = as.data.frame(entrez_id)
gene_df$gene_name <- row.names(gene_df)

doids <- select(
  x = HDO.db,
  keys = gene_df$entrez_id, keytype = "gene",
  columns = c("doid"))

colnames(doids) <- c("doid", "entrez_id")
full <- merge(doids, gene_df, on="entrez_id", all=TRUE)
full %>%
    group_by(gene_name) %>%
    summarize(n_doid = n()) -> n_doid

doid_degree <- merge(degree_df, n_doid, by="gene_name")
write.table(doid_degree,annotated_degree,sep="\t", row.names = FALSE)

