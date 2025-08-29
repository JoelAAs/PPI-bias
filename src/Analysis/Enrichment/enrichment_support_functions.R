library(org.Hs.eg.db)
library(clusterProfiler)


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
enrichment_DO_results <- function(genes){
  enrichTest <- EnrichDO::doEnrich(genes, method = "BH")
  # Multivalue cell values
  results <- as_tibble(enrichTest@enrich)
  #for (column_name in colnames(results)){
  #  if (type(results[, column_name]) == "list") {
  #    results[, column_name] = sapply(
  #        results[, column_name], function(x) paste(unlist(x), collapse=";")
  #        )
  #  }
  #}
  return(results)
}
enrichment_GO_results <- function(genes, ont="MF"){
  ego <- enrichGO(gene          = genes,
                  OrgDb         = org.Hs.eg.db,
                  ont           = ont, # Cellular components
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  return(ego@result)
}

get_enrichments <- function(genes) {
    go_enrichment <- enrichment_GO_results(genes)
    do_enrichment <- enrichment_DO_results(genes)
    return(
        list(
            go_results = go_enrichment,
            do_results = do_enrichment
            )
        )
}





