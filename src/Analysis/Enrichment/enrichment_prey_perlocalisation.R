home="/home/joel/Projects/cellular_compartment/"
source(paste0(home,"src/enrichment_support_functions.R"))
library(tidyverse)
library(ggVennDiagram)
library(magrittr)


### Functions
enrichment_per_bait_group <- function(df){
  i = 1
  go_list <- list()
  do_list <- list()
  for (localisation in unique(df$localisation_bait)) {
    bait_go_df <- df |> 
      filter(localisation_bait == localisation)
    go_result <- enrichment_GO_results(bait_go_df$entrez_id)
    go_result$localisation   <- localisation
    go_result$bait_N <- nrow(bait_go_df) 
    do_result <- enrichment_DO_results(bait_go_df$entrez_id)
    do_result$localisation   <- localisation
    do_result$bait_N <- nrow(bait_go_df)
    go_list[[i]] = go_result
    do_list[[i]] = do_result
      i = i+1
  }
  return(list(
    go_results = go_list, 
    do_results = do_list)
    )
}

# Run
prey_probabilies <- read.table(paste0(home, "work_folder/prey_probability/bait_all.csv"), header=T)
prey_probabilies |>
  filter(in_bioid=="True") |>
  filter(in_permutation > 0) -> prey_probabilies

prey_probabilies <- map_entrez(prey_probabilies, "gene_name_prey")
prey_probabilies$z <- (prey_probabilies$observed_value - prey_probabilies$probability_mean)/prey_probabilies$probability_std

ggplot(prey_probabilies,
       aes(x=z)) +
  geom_histogram() +
  facet_wrap(localisation_bait ~.)

prey_probabilies |>
  top_frac(0.1, z) -> top_5

prey_probabilies |>
  top_frac(0.1, -z) -> bot_5

top_enrichments <- enrichment_per_bait_group(top_5)
bot_enrichments = enrichment_per_bait_group(bot_5)

top_go_df <- bind_rows(top_enrichments$go_results)
top_do_df <- bind_rows(top_enrichments$do_results)

bot_go_df <- bind_rows(bot_enrichments$go_results)
bot_do_df <- bind_rows(bot_enrichments$do_results)

top_go_df %>% 
  filter(qvalue < 0.05) -> significant_enrichment_top
significant_enrichment_top %>% 
  {table(.$ID)} %>% 
  stack() -> significant_enrichment_top_count
colnames(significant_enrichment_top_count) <- c(
  "top_count",
  "go_term"
)

bot_go_df %>% 
  filter(qvalue < 0.05) -> significant_enrichment_bot
significant_enrichment_bot %>%  
  {table(.$ID)} %>% 
  stack() -> significant_enrichment_bot_count
colnames(significant_enrichment_bot_count) <- c(
  "bot_count",
  "go_term"
)

venn_enrich <- list(Top = significant_enrichment_top$ID, Bot = significant_enrichment_bot$ID)

go_vd <- ggVennDiagram(venn_enrich) +
  ggtitle("Significantly enriched CC GO terms camping top/bot 10% z-scores")

write.table(significant_enrichment_top, paste0(home,"work_folder/enrichment_analysis/enrichment/top_GO_z_enrichment.csv"), sep ="\t")
write.table(significant_enrichment_bot, paste0(home,"work_folder/enrichment_analysis/enrichment/bot_GO_z_enrichment.csv"), sep ="\t")