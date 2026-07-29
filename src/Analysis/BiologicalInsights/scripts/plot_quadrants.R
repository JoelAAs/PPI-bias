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
    log_deg_pos = log10(deg_pos + 1),
    log_deg_neg = log10(deg_neg + 1)
  )

g <- ggplot(df, aes(x = log_deg_neg, y = log_deg_pos)) +
  geom_hex(bins = 40) +
  facet_wrap(~ dataset_label, nrow = 1) +
  scale_fill_viridis_c(trans = "log10") +
  labs(
    title = "Degree quadrants",
    x = expression(log[10]("negative degree" + 1)),
    y = expression(log[10]("positive degree" + 1)),
    fill = "N proteins"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 9))

ggsave(snakemake@output[[1]], g, dpi = 300, height = 4, width = 9)
