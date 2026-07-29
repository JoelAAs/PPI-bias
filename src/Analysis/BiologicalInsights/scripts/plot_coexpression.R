library(ggplot2)
library(dplyr)
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

read_one <- function(path) {
  dataset <- sub(paste0("_", network_type, "_pairs\\.tsv$"), "", basename(path))
  read.table(path, sep = "\t", header = TRUE) %>%
    mutate(dataset = dataset)
}

df <- bind_rows(lapply(snakemake@input$pairs, read_one)) %>%
  mutate(
    dataset_label = factor(dataset_label(dataset), levels = c("Y2H", "MS", "Combined")),
    edge_set_label = case_when(
      edge_set == "hri" ~ "HRI",
      edge_set == "hrni" ~ "HRNI",
      edge_set == "random_universe_pairs" ~ "Random pairs",
      TRUE ~ edge_set
    )
  )

random_medians <- df %>%
  filter(edge_set == "random_universe_pairs") %>%
  group_by(dataset_label) %>%
  summarise(median_r = median(coexpr), .groups = "drop")

g <- ggplot(
  df %>% filter(edge_set != "random_universe_pairs"),
  aes(x = coexpr, fill = edge_set_label, color = edge_set_label)
) +
  geom_density(alpha = 0.4) +
  geom_vline(
    data = random_medians, aes(xintercept = median_r),
    linetype = "dashed", color = "grey40"
  ) +
  facet_wrap(~ dataset_label, nrow = 1) +
  scale_fill_manual(values = c("HRI" = "steelblue", "HRNI" = "darkorange")) +
  scale_color_manual(values = c("HRI" = "steelblue", "HRNI" = "darkorange")) +
  labs(
    title = "Co-expression by edge set (dashed line: median of random universe pairs)",
    x = "Pearson co-expression",
    y = "Density",
    fill = "Edge set",
    color = "Edge set"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 9), legend.position = "bottom")

ggsave(snakemake@output$density, g, dpi = 300, height = 4, width = 8)

## ---- sum_tests vs mean co-expression (per protein) ----
## Confound check: are heavily-tested proteins (protein_degrees.py's sum_tests, same
## edge population used here) more or less co-expressed with their partners on
## average, simply because they were tested more?
read_degree <- function(path) {
  dataset <- sub(paste0("_", network_type, "_degree\\.pq$"), "", basename(path))
  read_parquet(path) %>%
    select(uniprot_id, sum_tests) %>%
    mutate(dataset = dataset)
}

read_protein_coexpr <- function(path) {
  dataset <- sub(paste0("_", network_type, "_protein_mean_coexpr\\.tsv$"), "", basename(path))
  read.table(path, sep = "\t", header = TRUE) %>%
    mutate(dataset = dataset)
}

degree_df <- bind_rows(lapply(snakemake@input$degree, read_degree))
protein_coexpr_df <- bind_rows(lapply(snakemake@input$protein_coexpr, read_protein_coexpr))

sum_tests_df <- inner_join(degree_df, protein_coexpr_df, by = c("uniprot_id", "dataset")) %>%
  mutate(dataset_label = factor(dataset_label(dataset), levels = c("Y2H", "MS", "Combined")))

g_sum_tests <- ggplot(sum_tests_df, aes(x = sum_tests, y = mean_coexpr)) +
  geom_hex(bins = 50) +
  scale_x_log10() +
  scale_fill_continuous(type = "viridis") +
  facet_wrap(~ dataset_label, nrow = 1) +
  labs(
    title = "Per-protein test intensity (sum_tests) vs mean co-expression with tested partners",
    x = expression(log[10]("sum_tests")),
    y = "Mean Pearson co-expression (r) with tested partners",
    fill = "Count"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 9))

ggsave(snakemake@output$sum_tests_vs_coexpression, g_sum_tests, dpi = 300, height = 4, width = 9)
