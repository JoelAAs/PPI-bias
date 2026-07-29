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
  dataset <- sub(paste0("_", network_type, "_bins\\.tsv$"), "", basename(path))
  read.table(path, sep = "\t", header = TRUE) %>%
    mutate(dataset = dataset)
}

df <- bind_rows(lapply(snakemake@input, read_one)) %>%
  mutate(
    dataset_label = factor(dataset_label(dataset), levels = c("Y2H", "MS", "Combined")),
    degree_diff = deg_pos - deg_neg
  )

g <- ggplot(df, aes(x = n_pubmed + 1, y = degree_diff)) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_smooth(method = "loess", color = "steelblue", se = TRUE) +
  scale_x_log10() +
  facet_wrap(~ dataset_label, nrow = 1) +
  labs(
    title = "Degree difference vs research attention",
    x = expression(log[10]("n_pubmed" + 1)),
    y = "deg_pos - deg_neg"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 9))

ggsave(snakemake@output$diff_vs_references, g, dpi = 300, height = 4, width = 9)

g_hex <- ggplot(df, aes(x = deg_pos / sum_tests, y = deg_neg / sum_tests)) +
  geom_hex(bins = 50) +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_continuous(type = "viridis") +
  facet_wrap(~ dataset_label, nrow = 1) +
  labs(
    title = "Positive vs negative degree, normalized by test count",
    x = expression(log[10]("deg_pos / sum_tests")),
    y = expression(log[10]("deg_neg / sum_tests")),
    fill = "Count"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 9))

ggsave(snakemake@output$degree_vs_tests, g_hex, dpi = 300, height = 4, width = 9)

g_density <- ggplot(df, aes(x = sum_tests + 1)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  scale_x_log10() +
  facet_wrap(~ dataset_label, nrow = 1) +
  labs(
    title = "Test intensity (sum_tests) distribution",
    x = expression(log[10]("sum_tests" + 1)),
    y = "Density"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 9))

ggsave(snakemake@output$sum_tests_density, g_density, dpi = 300, height = 4, width = 9)


g_hex <- ggplot(df, aes(y = deg_pos + deg_neg, x = sum_tests)) +
  geom_hex(bins = 50) +
  scale_fill_continuous(type = "viridis") +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ dataset_label, nrow = 1) +
  labs(
    title = "Degree vs number of tests done on node",
    x = expression(log[10]("N tests")),
    fill = "Count",
    y = expression(log[10]("deg_neg + deg_pos")),
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 9))

ggsave(snakemake@output$test_vs_degree, g_hex, dpi = 300, height = 4, width = 9)