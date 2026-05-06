library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

protein_data <- read.csv(
  snakemake@input$per_protein_similarity,
  sep = "\t", header = TRUE
)
pair_data <- read.csv(
  snakemake@input$per_edge_similarity,
  sep = "\t", header = TRUE
)

wide <- protein_data %>%
  pivot_wider(
    names_from = set,
    values_from = bitscore_p_residue,
    names_glue = "{set}_bitscore_p_residue"
  )

wide$delta_random_negative <- wide$random_negative_bitscore_p_residue - wide$positive_bitscore_p_residue
wide$delta_negative <- wide$negative_bitscore_p_residue - wide$positive_bitscore_p_residue


wide_delta <- wide %>%
   filter(delta_negative != 0 | delta_random_negative != 0)
# wide_random <- wide %>%
#   filter(random_negative_bitscore_p_residue != 0 | positive_bitscore_p_residue != 0)

g_delta <- ggplot(
  wide_delta,
  aes(x = delta_random_negative, y = delta_negative)
) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  lims(x = c(-20, 20), y = c(-20, 20)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    title = "Per-protein sequence similarity delta to positives",
    x = "Δ avg bitscore per residue (random negative − positive)",
    y = "Δ avg bitscore per residue (negative − positive)",
    color = "Density"
  ) +
  theme_bw()

# g_rand <- ggplot(
#   wide_random,
#   aes(x = random_negative_bitscore_p_residue, y = positive_bitscore_p_residue)
# ) +
#   geom_hex(bins = 70) +
#   scale_fill_continuous(type = "viridis") +
#   lims(x = c(0, 10), y = c(0, 10)) +
#   geom_abline(intercept = 0, slope = 1, lintype = "dashed") +
#   labs(
#     title = "Non-zero per-protein similarity: positive vs random negative",
#     x = "Avg bitscore per protein (random negative edges)",
#     y = "Avg bitscore per protein (positive edges)",
#     color = "Density"
#   ) +
#   theme_bw()


# g_multi <- grid.arrange(g_neg, g_rand, ncol = 2)
ggsave(
  snakemake@output$avg_protein_similarity, g_delta,
  width = 6, height = 5, dpi = 150
)

g_pair <- ggplot(pair_data, aes(x = log(bitscore_p_residue), fill = set)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~set, ncol = 1) +
  labs(
    title = "Pair-level sequence similarity by set",
    x = "Non-zero Bitscore per residue",
    y = "Density",
    fill = "Set"
  ) +
  theme_bw() +
  theme(legend.position = "none")
ggsave(snakemake@output$edge_similarity_density, g_pair, width = 6, height = 7, dpi = 150)
