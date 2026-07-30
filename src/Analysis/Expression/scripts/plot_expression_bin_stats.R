library(ggplot2)
library(dplyr)
library(ggsignif)

## ---- Plot 1: co-expression by edge set (HRI vs HRNI), Wilcoxon bracket ----
df <- read.table(snakemake@input$hrni_vs_hri, sep = "\t", header = TRUE) %>%
  mutate(pos_label = factor(pos, levels = c(0, 1), labels = c("HRNI", "HRI")))

g_box <- ggplot(df, aes(x = pos_label, y = coexpr, fill = pos_label)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  geom_signif(
    comparisons = list(c("HRNI", "HRI")), test = wilcox.test,
    map_signif_level = FALSE, textsize = 3.5
  ) +
  scale_fill_manual(values = c("HRNI" = "darkorange", "HRI" = "steelblue")) +
  labs(
    title = "Co-expression by edge set",
    x = NULL, y = "Pearson co-expression"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10), legend.position = "none")

ggsave(snakemake@output$box_plot, g_box, dpi = 300, height = 4, width = 5)

## ---- Plot 2: odds ratio of expression-bin membership, HRI vs HRNI ----
or_df <- read.table(snakemake@input$expression_OR, sep = "\t", header = TRUE) %>%
  mutate(
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

g_bar <- ggplot(or_df, aes(x = class_label, y = log_or, fill = class_label)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(aes(y = text_y, label = sig), size = 5) +
  scale_fill_manual(values = c("Low" = "darkorange", "High" = "steelblue")) +
  labs(
    title = "Odds ratio of expression-bin membership: HRI vs HRNI",
    x = "Expression bin", y = expression(log("Odds ratio"))
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10), legend.position = "none")

ggsave(snakemake@output$or_barplot, g_bar, dpi = 300, height = 4, width = 5)
