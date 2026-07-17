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
  dataset <- sub(paste0("_", network_type, "\\.tsv$"), "", basename(path))
  read.table(path, sep = "\t", header = TRUE) %>%
    mutate(dataset = dataset)
}

df <- bind_rows(lapply(snakemake@input, read_one)) %>%
  mutate(
    dataset_label = factor(
      dataset_label(dataset),
      levels = c("Y2H", "MS", "Combined")
    ),
    log10_or = log10(odds_ratio),
    log10_ci_lo = log10(ci_lo),
    log10_ci_hi = log10(ci_hi)
  )

annotation_order <- df %>%
  group_by(annotation_type, annotation) %>%
  summarise(mean_log10_or = mean(log10_or), .groups = "drop") %>%
  arrange(annotation_type, mean_log10_or) %>%
  pull(annotation)

df <- df %>%
  mutate(annotation = factor(annotation, levels = annotation_order))

max_log10_ci_hi <- max(df$log10_ci_hi, na.rm = TRUE)

g <- ggplot(
  df,
  aes(x = log10_or, y = annotation, color = annotation_type)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = log10_ci_lo, xmax = log10_ci_hi), height = 0) +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(-max_log10_ci_hi, max_log10_ci_hi)) +
  scale_color_manual(
    values = c("GO" = "blue", "localisation" = "darkorange")
  ) +
  facet_wrap(~ dataset_label, nrow = 1) +
  labs(
    title = "Shared-annotation odds ratio",
    x = expression(log[10](odds~ratio)),
    y = "",
    color = "Annotation type"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 6),
    legend.position = "bottom"
  )

ggsave(snakemake@output[[1]], g, dpi = 300, height = 5, width = 5)
