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

dataset_colors <- c("Combined" = "blue", "MS" = "darkgreen", "Y2H" = "darkorange")

read_one <- function(path) {
  dataset <- sub(paste0("_", network_type, "_coessentiality_OR\\.tsv$"), "", basename(path))
  read.table(path, sep = "\t", header = TRUE) %>%
    mutate(dataset = dataset)
}

df <- bind_rows(lapply(snakemake@input$coessentiality_OR, read_one)) %>%
  mutate(
    dataset_label = factor(dataset_label(dataset), levels = c("Combined", "MS", "Y2H")),
    log10_or = log10(odds_ratio),
    log10_ci_lo = log10(ci_lo),
    log10_ci_hi = log10(ci_hi)
  )

g <- ggplot(df, aes(x = dataset_label, y = log10_or, fill = dataset_label)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = log10_ci_lo, ymax = log10_ci_hi), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = dataset_colors) +
  labs(
    title = "Co-essentiality: HRI vs HRNI",
    x = "Dataset", y = expression(log[10]("Odds ratio")), fill = "Dataset"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(snakemake@output[[1]], g, dpi = 300, height = 4, width = 3)
