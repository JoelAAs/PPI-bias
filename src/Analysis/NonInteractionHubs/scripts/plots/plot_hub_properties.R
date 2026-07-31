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
  ggsave(snakemake@output[[2]], g, dpi = 300, height = 3, width = 6)
  quit(save = "no")
}

df <- df %>%
  mutate(
    dataset_label = factor(dataset_label(dataset), levels = c("Y2H", "MS", "Combined")),
    sig_label = ifelse(!is.na(q_val) & q_val < 0.05, "*", "")
  )

dataset_colors <- c("Y2H" = "#2a78d6", "MS" = "#eb6834", "Combined" = "#1baf7a")

effect_plot <- function(sub_df, x_var, x_label, title) {
  dodge <- position_dodge(width = 0.6)
  ggplot(sub_df, aes(x = .data[[x_var]], y = feature, color = dataset_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(size = 2, position = dodge) +
    geom_text(
      aes(label = sig_label),
      position = dodge, hjust = -0.6, vjust = 0.35, size = 5, fontface = "bold",
      show.legend = FALSE
    ) +
    scale_color_manual(values = dataset_colors, name = "Dataset") +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.15))) +
    labs(title = title, subtitle = "* q < 0.05 (BH-corrected)", x = x_label, y = "") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 9),
      plot.subtitle = element_text(size = 8, color = "grey30")
    )
}

binary_df <- df %>% filter(feature_type == "binary") %>% mutate(log_or = log10(effect))
continuous_df <- df %>% filter(feature_type == "continuous")

g_binary <- effect_plot(binary_df, "log_or", "log10(Odds ratio)",
                         "Non-interaction hub vs actual-interactor reference: binary features")
g_continuous <- effect_plot(continuous_df, "effect", "% difference in mean",
                             "Non-interaction hub vs actual-interactor reference: continuous features")

ggsave(snakemake@output[[1]], g_binary, dpi = 300, height = 5, width = 7)
ggsave(snakemake@output[[2]], g_continuous, dpi = 300, height = 5, width = 7)
