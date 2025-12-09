library(HDO.db)
library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)
library(reshape2)

map_entrez <- function(gene_name_df, gene_column="gene_name") {
  entrez_ids <- mapIds(
    org.Hs.eg.db, keys = gene_name_df[, gene_column],
    column = "ENTREZID", keytype = "SYMBOL")
  stacked = stack(entrez_ids)
  stacked_df = as.data.frame(stacked)
  colnames(stacked_df) <- c("entrez_id", gene_column)
  gene_name_df$entrez_id = stacked_df$entrez_id
  return(gene_name_df)
}

dplot_go <- function(df) {
  df$GeneRatio <- sapply(df$GeneRatio, function(x) eval(parse(text=x)))
  df$significant <- df$qvalue < 0.05
  df <- df %>%
    mutate(Description = fct_reorder(Description, GeneRatio))
  p <- ggplot(
    df,
    aes(y=Description,
        x=dataset,
        color=-log10(qvalue),
        size=GeneRatio,
        shape=significant)
    ) + geom_point()
  return(p)
}

dplot_do <- function(df, color="significant") {
  df$significant <- df$qvalue < 0.05
  df$GeneRatio <- sapply(df$GeneRatio, function(x) eval(parse(text=x)))
  df <- df %>%
    mutate(Description = fct_reorder(Description, GeneRatio))
  p <- ggplot(
    df,
    aes(y=Description,
        x=dataset,
        color=!!sym(color),
        size=GeneRatio)
  ) + geom_point()
  return(p)
}


get_enrich_go  <- function(df, n, val_column="norm_degree") {
  df %>% 
    slice_max(!!sym(val_column), n=n) %>%
    pull(entrez_id) -> genes
  
    enrichGO(
      gene          = genes,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP", 
      pAdjustMethod = "BH",
      qvalueCutoff  = 1,
      readable      = TRUE)@result %>%
    return()
}

get_enrich_do  <- function(df, n, val_column="norm_degree") {
  df %>% 
    slice_max(!!sym(val_column), n=n) %>%
    pull(entrez_id) -> genes
  clusterProfiler::enrichDO(genes, qvalueCutoff = 1, pvalueCutoff = 1)@result %>%
    return()
}


get_top_id <- function(enrich_df, n, col="ID", value = "qvalue", significant=T) {
  if (significant) {
    enrich_df %>% filter(!!sym(value) < 0.05) -> enrich_df
  }
  enrich_df %>% 
  slice_min(!!sym(value), n=n) %>%
  pull(!!sym(col)) %>% return()
}

get_gene_doid_overlap <- function(enrich_results, top_genes){
  enrich_results@result$geneID_split <- sapply(enrich_results@result$geneID, function(x) str_split(x, "/"))
  i = 1
  n_doids = rep(0, length(enrich_results@gene))
  for (gene in enrich_results@gene) {
    n_doids[i] = sum(sapply(enrich_results$geneID_split, function(x) gene %in% x))
    i = i + 1
  }
  data.frame(
    entrez_id = enrich_results@gene,
    n_doid = n_doids,
    among_top = sapply(enrich_results@gene, function(x) x %in% top_genes)
  ) %>% return()
}
  
get_norm_degree_delta <- function(df1, df2, merged_df, target_col, col_name, missing = T){
  cols <- c("gene_name", target_col)
  delta = merge(df1[, cols], df2[, cols], by="gene_name", all.x=TRUE, suffixes=c("_1", "_2"))
  if (missing){
    delta[is.na(delta[, paste0(target_col, "_2")]), paste0(target_col, "_2")] <- 0
  }
  delta[, col_name] <- delta[, paste0(target_col, "_1")] - delta[, paste0(target_col, "_2")]
  merged_df = merge(merged_df, delta[, c("gene_name", col_name)])
  return(merged_df)
}

args <- commandArgs(trailingOnly = TRUE)

degree_file <- args[1]
go_enrichment_bait  <- args[2]
go_enrichment_prey  <- args[3]
do_enrichment_bait  <- args[4]
do_enrichment_prey  <- args[5]
df_degree = read.table(degree_file, sep ="\t", header=T) %>%
  map_entrez() %>% filter(!is.na(entrez_id))

n_genes = 50
## GO degree
df_degree %>% get_enrich_go(n_genes, "degree_bait") %>% filter(qvalue < 0.05) -> go_degree_bait
df_degree %>% get_enrich_go(n_genes, "degree_prey") %>% filter(qvalue < 0.05) -> go_degree_prey

## DO degree
df_degree %>% get_enrich_do(n_genes, "degree_bait") %>% filter(qvalue < 0.05) -> do_degree_bait
df_degree %>% get_enrich_do(n_genes, "degree_prey") %>% filter(qvalue < 0.05) -> do_degree_prey

write.table(go_degree_bait,go_enrichmentt_bait,sep="\t", row_name=F)
write.table(go_degree_prey,go_enrichment_prey,sep="\t", row_name=F)
write.table(do_degree_bait,do_enrichmentt_bait,sep="\t", row_name=F)
write.table(do_degree_prey,do_enrichment_prey,sep="\t", row_name=F)

