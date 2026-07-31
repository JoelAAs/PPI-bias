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
  dataset <- sub(paste0("_", network_type, "_degree\\.pq$"), "", basename(path))
  read.table(path, sep = "\t", header = TRUE) %>%
    mutate(dataset = dataset)
}

df <- bind_rows(lapply(snakemake@input, read_one)) %>%
  mutate(
    dataset_label = factor(dataset_label(dataset), levels = c("Y2H", "MS", "Combined")),
    log_deg_pos = log10(deg_pos),
    log_deg_neg = log10(deg_neg)
  )

g <- ggplot(df, aes(x = log_deg_neg, y = log_deg_pos)) +
  geom_hex(bins = 40) +
  facet_wrap(~ dataset_label, nrow = 1) +
  scale_fill_viridis_c(trans = "log10") +
  labs(
    title = "Degree/non-degree distribution",
    x = expression(log[10]("negative degree")),
    y = expression(log[10]("positive degree")),
    fill = "N proteins"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 9))

ggsave(snakemake@output[[1]], g, dpi = 300, height = 4, width = 9)


g <- ggplot(df, aes(x = n_pubmed + 1, y = deg_pos - deg_neg)) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_smooth(method = "loess", color = "steelblue", se = TRUE) +
  scale_x_log10() +
  facet_wrap(~ dataset_label, nrow = 1) +
  labs(
    title = "Degree difference vs research attention",
    x = expression(log[10]("n_pubmed" + 1)),
    y = "Interaction Degree - Non-interaction Degree"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 9))

ggsave(snakemake@output[[2]], g, dpi = 300, height = 4, width = 9)
