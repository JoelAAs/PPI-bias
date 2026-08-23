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

df <- df %>%
  mutate(
    dataset_label = factor(dataset_label(dataset), levels = c("Combined", "MS", "Y2H")),
    sig_label = ifelse(!is.na(q_val) & q_val < 0.05, "*", "")
  )

dataset_colors <- c("Combined" = "blue", "MS" = "darkgreen", "Y2H" = "darkorange")

effect_plot <- function(sub_df, x_var, x_label, title) {
  dodge <- position_dodge(width = 0.6)
  ggplot(sub_df, aes(x = .data[[x_var]], y = feature, fill = dataset_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_col(width = 0.6, position = dodge) +
    geom_text(
      aes(label = sig_label, hjust = ifelse(.data[[x_var]] >= 0, -0.6, 1.6)),
      position = dodge, vjust = 0.35,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = dataset_colors, name = "Dataset") +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.15))) +
    labs(title = title, x = x_label, y = "") +
    theme_bw() +
    theme(legend.position = "bottom")
}

binary_df <- df %>% filter(feature_type == "binary") %>% mutate(log_or = log10(effect))
continuous_df <- df %>% filter(feature_type == "continuous")

g_binary <- effect_plot(binary_df, "log_or", expression(log[10]("Odds ratio")),
                         "NI hub vs I hub: binary features")
g_continuous <- effect_plot(continuous_df, "effect", "% difference in mean",
                             "NI hub vs I hub: continuous features")

ggsave(snakemake@output[[1]], g_binary, height = 4, width = 5)
ggsave(snakemake@output[[2]], g_continuous, height = 4, width = 5)
