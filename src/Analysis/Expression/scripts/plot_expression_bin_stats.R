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

dataset_colors <- c("Combined" = "blue", "MS" = "darkgreen", "Y2H" = "darkorange")

read_hri_hrni <- function(path) {
  dataset <- sub(paste0("_", network_type, "_hri_hnri\\.tsv$"), "", basename(path))
  read.table(path, sep = "\t", header = TRUE) %>%
    mutate(dataset = dataset)
}

group_levels <- unlist(lapply(c("Combined", "MS", "Y2H"), function(d) paste(d, c("HRNI", "HRI"), sep = "\n")))

df <- bind_rows(lapply(snakemake@input$hrni_vs_hri, read_hri_hrni)) %>%
  mutate(
    dataset_label = factor(dataset_label(dataset), levels = c("Combined", "MS", "Y2H")),
    pos_label = factor(pos, levels = c(0, 1), labels = c("HRNI", "HRI")),
    group_label = factor(paste(dataset_label, pos_label, sep = "\n"), levels = group_levels)
  )

g_box <- ggplot(df, aes(x = group_label, y = coexpr, fill = dataset_label)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  geom_signif(
    comparisons = list(
      c("Combined\nHRNI", "Combined\nHRI"),
      c("MS\nHRNI", "MS\nHRI"),
      c("Y2H\nHRNI", "Y2H\nHRI")
    ),
    test = wilcox.test, map_signif_level = FALSE, textsize = 3.5
  ) +
  scale_fill_manual(values = dataset_colors) +
  labs(
    title = "Co-expression, HRI/HRNI sets",
    x = NULL, y = "Pearson co-expression", fill = "Dataset"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10), legend.position = "bottom")

ggsave(snakemake@output$box_plot, g_box, dpi = 300, height = 4, width = 8)

read_expression_or <- function(path) {
  dataset <- sub(paste0("_", network_type, "_expression_OR\\.tsv$"), "", basename(path))
  read.table(path, sep = "\t", header = TRUE) %>%
    mutate(dataset = dataset)
}

or_df <- bind_rows(lapply(snakemake@input$expression_OR, read_expression_or)) %>%
  mutate(
    dataset_label = factor(dataset_label(dataset), levels = c("Combined", "MS", "Y2H")),
    log_or = log10(OR),
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    text_y = log_or + ifelse(log_or >= 0, 0.05, -0.05),
    class_label = factor(class, levels = c("low", "high"), labels = c("Bottom 15%", "Top 15%"))
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
    title = "Expression-bin: HRI vs HRNI",
    x = "Expression bin", y = expression(log[10]("Odds ratio")), fill = "Dataset"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10), legend.position = "bottom")

ggsave(snakemake@output$or_barplot, g_bar, dpi = 300, height = 4, width = 4)
