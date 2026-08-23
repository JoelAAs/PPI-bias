library(ggplot2)
library(dplyr)
library(tidyr)

network_type <- snakemake@wildcards$network_type

dataset_label <- function(dataset) {
  case_when(
    dataset == "ms" ~ "MS",
    dataset == "y2h" ~ "Y2H",
    TRUE ~ dataset
  )
}

dataset_colors <- c("MS" = "darkgreen", "Y2H" = "darkorange")

read_one <- function(path) {
  dataset <- sub(paste0("_", network_type, "_detection\\.csv$"), "", basename(path))
  read.table(path, sep = "\t", header = TRUE) %>%
    mutate(dataset = dataset)
}

df <- bind_rows(lapply(snakemake@input, read_one)) %>%
  mutate(
    dataset_label = factor(dataset_label(dataset), levels = c("MS", "Y2H")),
    Category = factor(Category, levels = c("large (strong)", "medium", "small (weak)", "Unknown"))
  )

df_plot <- df %>%
  select(Category, dataset_label, pair_detection_ratio, global_detection_ratio) %>%
  pivot_longer(
    cols = c(pair_detection_ratio, global_detection_ratio),
    names_to = "metric", values_to = "value"
  )

g <- ggplot(
  df_plot,
  aes(x = Category, y = value, color = dataset_label,
      group = interaction(dataset_label, metric), shape = metric)
) +
  geom_point(size = 2) +
  geom_line() +
  theme_bw() +
  labs(
    title = "Detection ratio by interface size: pair-wise vs global",
    x = "Interface size",
    y = "Detection ratio",
    shape = "Metric",
    color = "Dataset"
  ) +
  scale_color_manual(values = dataset_colors) +
  scale_shape_manual(
    values = c(pair_detection_ratio = 15, global_detection_ratio = 8),
    labels = c(pair_detection_ratio = "Pair-wise", global_detection_ratio = "Pooled")
  ) +
  theme(legend.position = "bottom")

ggsave(snakemake@output[[1]], g, dpi = 300, height = 4, width = 6)
