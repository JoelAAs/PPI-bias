library(ggplot2)
library(dplyr)

bp_df <- read.table(snakemake@input$bp_df, sep = "\t", header = TRUE)

bp_df %>%
  filter(detection_method %in% snakemake@params$ms_methods) %>%
  mutate(
    uniprot_id_bait = sapply(
      uniprot_id_bait, function(x) strsplit(x, "-")[[1]][1]),
    uniprot_id_prey = sapply(
      uniprot_id_prey, function(x) strsplit(x, "-")[[1]][1])
  ) %>%
  filter(uniprot_id_bait != uniprot_id_prey) -> bp_df

study_counts <- bp_df %>%
  mutate(study = paste(pubmed_id, detection_method, sep = "_")) %>%
  group_by(study) %>%
  summarise(
    n_bait = n_distinct(uniprot_id_bait),
    n_prey = n_distinct(uniprot_id_prey),
    .groups = "drop"
  )

g <- ggplot(study_counts, aes(x = n_bait, y = n_prey)) +
  geom_point(alpha = 0.4) +
  geom_smooth(
    method = "nls",
    formula = y ~ SSmicmen(x, Vm, K),  # Michaelis-Menten: y = Vm * x / (K + x)
    se = FALSE,
    color = "steelblue"
  ) +
  labs(x = "Unique baits", y = "Unique prey") +
  theme_bw()

ggsave(snakemake@output$plot, g, dpi = 300, height = 5, width = 6)
