library(HDO.db)
library(tidyverse)
library(org.Hs.eg.db)

### args
args <- commandArgs(trailingOnly = TRUE)

input_pod  <- args[1]
output_pod <- args[2]

### Input
pod_df = read.table(
  input_pod,
  sep="\t",
  header=T
)

### Processing
all_gene <- unique(c(pod_df$gene_name_bait, pod_df$gene_name_prey))
entrez_id <- mapIds(
  org.Hs.eg.db, keys = all_gene,
  column = "ENTREZID", keytype = "SYMBOL")

gene_df = as.data.frame(entrez_id)
gene_df$gene_name <- row.names(gene_df)

doids <- select(
  x = HDO.db,
  keys = gene_df$entrez_id, keytype = "gene",
  columns = c("doid"))

colnames(doids) <- c("doid", "entrez_id")
full <- merge(doids, gene_df, on="entrez_id", all=TRUE)

write.table(full,output_pod,sep="\t", row.names = FALSE)
