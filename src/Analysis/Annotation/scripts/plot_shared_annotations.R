library(ggplot2)
library(dplyr)

df <- read.table(snakemake@input[[1]], sep = "\t", header = TRUE)

dataset_label <- case_when(
  snakemake@wildcards$dataset == "flat" ~ "Combined",
  snakemake@wildcards$dataset == "ms" ~ "MS",
  snakemake@wildcards$dataset == "y2h" ~ "Y2H",
  TRUE ~ snakemake@wildcards$dataset
)

df <- df %>%
#  filter(!(ci_lo <= 1 & ci_hi >= 1)) %>%
  arrange(annotation_type, odds_ratio) %>%
  mutate(annotation = factor(annotation, levels = annotation))

g <- ggplot(
  df,
  aes(x = odds_ratio, y = annotation, color = annotation_type)
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_color_manual(
    values = c("GO" = "blue", "localisation" = "darkorange")
  ) +
  labs(
    title = paste(
      "Shared-annotation odds ratio:", dataset_label
    ),
    x = "Odds ratio (log scale)",
    y = "Annotation",
    color = "Annotation type"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )

ggsave(snakemake@output[[1]], g, dpi = 300, height = 4, width = 4)
