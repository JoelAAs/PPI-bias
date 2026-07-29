# Module 3 - GO over-representation per degree quadrant (BiologicalInsights plan
# §8.3), with the TESTED PROTEIN UNIVERSE as background (not the whole proteome -
# that would produce generic terms reflecting only which proteins get screened).
# go_enrichment_source=config["clusterprofiler"] here: gprofiler2 is not installed in
# the do_enrichment env (only clusterProfiler, org.Hs.eg.db, HDO.db, GOfuncR are, per
# envs/setup_R_do_enrichment.R) so clusterProfiler::enrichGO is used, keyType="UNIPROT"
# (org.Hs.eg.db supports this directly - no extra ID mapping step needed).
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

min_gs_size <- 10
min_quadrant_size <- 5

bins_df <- read.table(snakemake@input[[1]], sep = "\t", header = TRUE)
universe <- unique(bins_df$uniprot_id)
cat(sprintf("Background universe: %d tested proteins\n", length(universe)),
    file = stderr())

quadrants <- unique(bins_df$quadrant)
results <- list()
for (q in quadrants) {
  genes <- bins_df$uniprot_id[bins_df$quadrant == q]
  if (length(genes) < min_quadrant_size) {
    cat(sprintf("Skipping quadrant %s: only %d proteins (< %d)\n",
                q, length(genes), min_quadrant_size), file = stderr())
    next
  }
  ego <- tryCatch(
    enrichGO(
      gene = genes,
      universe = universe,
      OrgDb = org.Hs.eg.db,
      keyType = "UNIPROT",
      ont = "BP",
      pAdjustMethod = "BH",
      minGSSize = min_gs_size
    ),
    error = function(e) {
      cat(sprintf("enrichGO failed for quadrant %s: %s\n", q, conditionMessage(e)),
          file = stderr())
      NULL
    }
  )
  if (!is.null(ego) && nrow(ego@result) > 0) {
    res <- as.data.frame(ego@result)
    res$quadrant <- q
    res$n_genes_in_quadrant <- length(genes)
    results[[q]] <- res
  }
}

if (length(results) > 0) {
  final <- bind_rows(results)
} else {
  final <- data.frame(
    quadrant = character(), n_genes_in_quadrant = integer(),
    ID = character(), Description = character(), GeneRatio = character(),
    BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
    qvalue = numeric(), geneID = character(), Count = integer()
  )
}
write.table(final, snakemake@output[[1]], sep = "\t", row.names = FALSE, quote = FALSE)
