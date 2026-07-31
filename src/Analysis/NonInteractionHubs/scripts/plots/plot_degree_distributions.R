library(ggplot2)
library(dplyr)
library(tidyr)
library(arrow)

network_type <- snakemake@wildcards$network_type

dataset_label <- function(dataset) {
  case_when(
    dataset == "flat" ~ "Combined",
    dataset == "ms" ~ "MS",
    dataset == "y2h" ~ "Y2H",
    TRUE ~ dataset
  )
}

extract_dataset <- function(path, suffix) {
  sub(paste0("_", network_type, suffix, "$"), "", basename(path))
}

degree_df <- bind_rows(lapply(snakemake@input$degree, function(path) {
  read_parquet(path) %>%
    mutate(dataset = extract_dataset(path, "_degree\\.pq"))
})) %>%
  mutate(dataset_label = factor(dataset_label(dataset), levels = c("Combined", "MS", "Y2H")))

bin_summary_df <- bind_rows(lapply(snakemake@input$bin_summary, function(path) {
  read.table(path, sep = "\t", header = TRUE) %>%
    mutate(dataset = extract_dataset(path, "_degree_bin_summary\\.tsv"))
})) %>%
  mutate(dataset_label = factor(dataset_label(dataset), levels = c("Combined", "MS", "Y2H")))

boundaries <- bin_summary_df %>%
  filter(bin %in% c("low", "medium")) %>%
  select(dataset_label, bin, max_deg_neg, max_deg_pos) %>%
  pivot_longer(
    cols = c(max_deg_neg, max_deg_pos),
    names_to = "degree_type", values_to = "boundary",
    names_prefix = "max_"
  )

long_degree <- degree_df %>%
  select(dataset_label, deg_neg, deg_pos) %>%
  pivot_longer(cols = c(deg_neg, deg_pos), names_to = "degree_type", values_to = "degree")

g_dist <- ggplot(long_degree, aes(x = degree + 1)) +
  geom_histogram(bins = 60) +
  geom_vline(
    data = boundaries,
    aes(xintercept = boundary + 1, linetype = bin),
    color = "firebrick"
  ) +
  scale_x_log10() +
  facet_grid(degree_type ~ dataset_label, scales = "free_y") +
  labs(
    title = "Degree distributions with low/medium/high bin limits",
    x = expression(log[10]("degree" + 1)),
    y = "N proteins",
    linetype = "Bin boundary"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 9))

ggsave(snakemake@output$distributions, g_dist, dpi = 300, height = 6, width = 9)

## ---- Plot 2: heatmap of deg_neg_bin x deg_pos_bin overlap ----
bin_levels <- c("low", "medium", "high")
crosstab <- degree_df %>%
  mutate(
    deg_neg_bin = factor(deg_neg_bin, levels = bin_levels),
    deg_pos_bin = factor(deg_pos_bin, levels = bin_levels)
  ) %>%
  count(dataset_label, deg_neg_bin, deg_pos_bin, .drop = FALSE)

g_heat <- ggplot(crosstab, aes(x = deg_neg_bin, y = deg_pos_bin, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), color = "white", size = 3) +
  scale_fill_viridis_c() +
  facet_wrap(~ dataset_label, nrow = 1) +
  labs(
    title = "Shared proteins between degree categories (deg_neg_bin x deg_pos_bin)",
    x = "deg_neg bin",
    y = "deg_pos bin",
    fill = "N proteins"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 9))

ggsave(snakemake@output$heatmap, g_heat, dpi = 300, height = 4, width = 9)
