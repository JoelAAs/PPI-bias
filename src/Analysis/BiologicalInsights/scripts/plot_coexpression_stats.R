library(ggplot2)
library(dplyr)

network_type <- snakemake@wildcards$network_type

dataset_label <- function(dataset) {
  case_when(
    dataset == "flat" ~ "Combined",
    dataset == "ms" ~ "MS",
    dataset == "y2h" ~ "Y2H",
    TRUE ~ dataset
  )
}

read_one <- function(path) {
  dataset <- sub(paste0("_", network_type, "_summary\\.tsv$"), "", basename(path))
  read.table(path, sep = "\t", header = TRUE) %>%
    mutate(dataset = dataset)
}

df <- bind_rows(lapply(snakemake@input$summary, read_one)) %>%
  mutate(dataset_label = factor(dataset_label(dataset), levels = c("Y2H", "MS", "Combined")))

## ---- Plot 1: odds ratio of "high" co-expression membership, HRI vs HRNI ----
or_df <- df %>% filter(analysis == "or_high_membership_hri_vs_hrni")

g_or <- ggplot(or_df, aes(x = dataset_label, y = effect)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), color = "steelblue", size = 0.7) +
  labs(
    title = "Odds ratio of \"high\" co-expression membership: HRI vs HRNI",
    subtitle = "Points: observed OR, bars: 95% cluster-bootstrap CI (dashed line: OR = 1, no association)",
    x = "Dataset",
    y = "Odds ratio"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

ggsave(snakemake@output$or_high, g_or, dpi = 300, height = 4, width = 6)

## ---- Plot 2: HRNI summed co-expression vs uniform-random null ----
sum_df <- df %>% filter(analysis == "summed_coexpression_vs_random")

g_sum <- ggplot(sum_df, aes(x = dataset_label)) +
  geom_pointrange(
    aes(y = null_mean, ymin = ci_lo, ymax = ci_hi, color = "Null (uniform random pairs)"),
    position = position_nudge(x = 0.15), size = 0.6
  ) +
  geom_point(
    aes(y = observed_sum, color = "Observed HRNI"),
    position = position_nudge(x = -0.15), size = 3
  ) +
  scale_color_manual(values = c("Observed HRNI" = "darkorange", "Null (uniform random pairs)" = "grey40")) +
  labs(
    title = "HRNI summed co-expression vs random pairs",
    subtitle = "Null range: 95th percentile interval over uniform random gene-pair draws (same n as HRNI)",
    x = "Dataset",
    y = "Summed co-expression (r)",
    color = NULL
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8), legend.position = "bottom")

ggsave(snakemake@output$summed, g_sum, dpi = 300, height = 4, width = 6)
