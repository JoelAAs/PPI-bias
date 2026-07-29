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
  dataset <- sub(paste0("_", network_type, "_enrichment\\.tsv$"), "", basename(path))
  read.table(path, sep = "\t", header = TRUE) %>%
    mutate(dataset = dataset)
}

df <- bind_rows(lapply(snakemake@input, read_one))

if (nrow(df) == 0 || all(is.na(df$effect))) {
  g <- ggplot() +
    annotate("text", x = 0, y = 0,
             label = "No hub-property results (hub set was empty upstream)", size = 4) +
    theme_void()
  ggsave(snakemake@output[[1]], g, dpi = 300, height = 3, width = 6)
  quit(save = "no")
}

df <- df %>%
  mutate(
    dataset_label = factor(dataset_label(dataset), levels = c("Y2H", "MS", "Combined")),
    log_effect = ifelse(feature_type == "binary", log10(effect), effect),
    log_ci_lo = ifelse(feature_type == "binary", log10(ci_lo), ci_lo),
    log_ci_hi = ifelse(feature_type == "binary", log10(ci_hi), ci_hi),
    sig_label = ifelse(!is.na(q_val) & q_val < 0.05, "*", "")
  )

g <- ggplot(df, aes(x = log_effect, y = feature, color = feature_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = log_ci_lo, xmax = log_ci_hi), height = 0) +
  geom_point(size = 2) +
  geom_text(
    aes(x = log_ci_hi, label = sig_label),
    hjust = -0.4, vjust = 0.35, size = 5, fontface = "bold", show.legend = FALSE
  ) +
  facet_wrap(~ dataset_label, nrow = 1) +
  scale_color_manual(values = c("binary" = "steelblue", "continuous" = "darkorange")) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(
    title = "No-interaction hub properties vs actual-interactor reference (binary: log10 OR; continuous: % diff)",
    subtitle = "* q < 0.05 (BH-corrected)",
    x = "Effect size",
    y = "",
    color = "Feature type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 9),
    plot.subtitle = element_text(size = 8, color = "grey30"),
    legend.position = "bottom"
  )

ggsave(snakemake@output[[1]], g, dpi = 300, height = 4, width = 8)
