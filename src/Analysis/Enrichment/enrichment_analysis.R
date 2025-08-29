source("src/enrichment_support_functions.R")
library(tidyverse)
library(ggVennDiagram)

#### functions
get_top_rows <- function(gene_name_df, n_top=100){
  entrez_id_df <- map_entrez(gene_name_df, "gene_name")
  gene_name_df |> 
    group_by(gene_name) |>
    summarise(n = n()) |>
    ungroup() -> df_bait_freq_df
  entrez_id_df <- entrez_id_df[!duplicated(entrez_id_df),]
  df_bait_freq_df <- merge(df_bait_freq_df, entrez_id_df, on="gene_name",how="inner")
  
  df_bait_freq_df$freq <- df_bait_freq_df$n/sum(df_bait_freq_df$n)
  df_bait_freq_df_top <- df_bait_freq_df |> top_n(n_top, freq)
  return(df_bait_freq_df_top)
}

#### args
args <- commandArgs(trailingOnly = TRUE)

bioid_bait_input                <- args[1]
ms_bait_input                   <- args[2]
bioid_bait_enrichment_go_output <- args[3]
bioid_bait_enrichment_do_output <- args[4]
ms_bait_enrichment_go_output    <- args[5]
ms_bait_enrichment_do_output    <- args[6]
venn_plot_go                    <- args[7]
venn_plot_doid                  <- args[8]
venn_plot_bait                  <- args[9]
bait_count                      <- as.numeric(args[10])

### run
df_bioid_baits <- read.table(bioid_bait_input, header=T)
df_ms_baits <- read.table(ms_bait_input, header=T)

# bioid
df_bioid_baits_top <- get_top_rows(df_bioid_baits, bait_count)
bioid_enrichments <- get_enrichments(df_bioid_baits_top$entrez_id)
write.table(bioid_enrichments$go_results, bioid_bait_enrichment_go_output, sep ="\t")
write.table(bioid_enrichments$do_results, bioid_bait_enrichment_do_output, sep ="\t")

# ms
df_ms_baits_top <- get_top_rows(df_ms_baits, bait_count)
ms_enrichments <- get_enrichments(df_ms_baits_top$entrez_id)
write.table(ms_enrichments$go_results, ms_bait_enrichment_go_output, sep ="\t")
write.table(ms_enrichments$do_results, ms_bait_enrichment_do_output, sep ="\t")

venn_go <- list(Bioid = bioid_enrichments$go_results$ID, MS = ms_enrichments$go_results$ID)
venn_doid <- list(Bioid = bioid_enrichments$do_results$DOID, MS = ms_enrichments$do_results$DOID)
venn_baits <- list(Bioid = df_bioid_baits_top$gene_name, MS = df_ms_baits_top$gene_name)

go_vd <- ggVennDiagram(venn_go) +
  ggtitle("Number of shared significantly enriched CC GO terms")
ggsave(venn_plot_go, go_vd, dpi=300)
doid_vd <- ggVennDiagram(venn_doid) +
  ggtitle("Number of shared significantly enriched DOIDs")
ggsave(venn_plot_doid, doid_vd, dpi=300)
bait_vd <- ggVennDiagram(venn_baits) +
  ggtitle(paste("Number of shared baits among the top", bait_count))
ggsave(venn_plot_bait, bait_vd, dpi=300)

