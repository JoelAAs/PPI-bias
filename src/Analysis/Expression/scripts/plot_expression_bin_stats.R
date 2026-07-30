library(ggplot2)
library(dplyr)
library(ggsignif)

network_type <- snakemake@wildcards$network_type

dataset_label <- function(dataset) {
  case_when(
    dataset == "flat" ~ "Combined",
    dataset == "ms" ~ "MS",
    dataset == "y2h" ~ "Y2H",
    TRUE ~ dataset
  )
}

dataset_colors <- c("Y2H" = "#2a78d6", "MS" = "#eb6834", "Combined" = "#1baf7a")

## ---- Plot 1: co-expression by edge set (HRI vs HRNI), Wilcoxon bracket ----
read_hri_hrni <- function(path) {
  dataset <- sub(paste0("_", network_type, "_hri_hnri\\.tsv$"), "", basename(path))
  read.table(path, sep = "\t", header = TRUE) %>%
    mutate(dataset = dataset)
}

df <- bind_rows(lapply(snakemake@input$hrni_vs_hri, read_hri_hrni)) %>%
  mutate(
    dataset_label = factor(dataset_label(dataset), levels = c("Y2H", "MS", "Combined")),
    pos_label = factor(pos, levels = c(0, 1), labels = c("HRNI", "HRI"))
  )

g_box <- ggplot(df, aes(x = pos_label, y = coexpr, fill = dataset_label)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7, position = position_dodge(width = 0.75)) +
  geom_signif(
    comparisons = list(c("HRNI", "HRI")), test = wilcox.test,
    map_signif_level = FALSE, textsize = 3.5
  ) +
  scale_fill_manual(values = dataset_colors) +
  labs(
    title = "Co-expression by edge set",
    x = NULL, y = "Pearson co-expression", fill = "Dataset"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10), legend.position = "bottom")

ggsave(snakemake@output$box_plot, g_box, dpi = 300, height = 4, width = 6)

## ---- Plot 2: odds ratio of expression-bin membership, HRI vs HRNI ----
read_expression_or <- function(path) {
  dataset <- sub(paste0("_", network_type, "_expression_OR\\.tsv$"), "", basename(path))
  read.table(path, sep = "\t", header = TRUE) %>%
    mutate(dataset = dataset)
}

or_df <- bind_rows(lapply(snakemake@input$expression_OR, read_expression_or)) %>%
  mutate(
    dataset_label = factor(dataset_label(dataset), levels = c("Y2H", "MS", "Combined")),
    log_or = log(OR),
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    text_y = log_or + ifelse(log_or >= 0, 0.05, -0.05),
    class_label = factor(class, levels = c("low", "high"), labels = c("Low", "High"))
  )

g_bar <- ggplot(or_df, aes(x = class_label, y = log_or, fill = dataset_label)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(
    aes(y = text_y, label = sig),
    position = position_dodge(width = 0.7), size = 5
  ) +
  scale_fill_manual(values = dataset_colors) +
  labs(
    title = "Odds ratio of expression-bin membership: HRI vs HRNI",
    x = "Expression bin", y = expression(log("Odds ratio")), fill = "Dataset"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10), legend.position = "bottom")

ggsave(snakemake@output$or_barplot, g_bar, dpi = 300, height = 4, width = 6)
