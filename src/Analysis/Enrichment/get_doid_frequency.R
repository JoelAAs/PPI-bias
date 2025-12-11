library(HDO.db)
library(tidyverse)
library(org.Hs.eg.db)

### args
args <- commandArgs(trailingOnly = TRUE)

input_degree  <- args[1]
n_top         <- as.numeric(args[2])
val_column        <- args[3]
doid_count_file <- args[4]
doid_annotation_file <- args[5]


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
degree_df <- merge(degree_df, gene_df, on="entrez_id",left_all=TRUE)
degree_df %>%
    slice_max(!!sym(val_column), n=n_top,with_ties = FALSE) -> degree_df
degree_df %>% pull(entrez_id) -> genes



doids <- select(
  x = HDO.db,
  keys = genes, keytype = "gene",
  columns = c("doid"))

colnames(doids) <- c("doid", "entrez_id")
degree_df_doid <- merge(degree_df, doids, on="entrez_id", left_all=TRUE)
degree_df_doid %>%
    group_by(doid) %>%
    summarize(
        doid_frequency = n()/n_top
    ) -> doid_count

write.table(doid_count,doid_count_file,sep="\t", row.names = FALSE)
write.table(degree_df_doid,doid_annotation_file,sep="\t", row.names = FALSE)

