library(ggplot2)

df <- read.table(snakemake@input$ms_sticky_proteins, sep = "\t", header = TRUE)

g <- ggplot(df, aes(x = dispersion_ratio)) +
  geom_histogram(bins = 50, fill = "forestgreen", color = "white") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  labs(
    x = "MS dispersion ratio: Var(observed) / Var(expected under uniform rate)",
    y = "Number of prey proteins"
  ) +
  theme_bw()

ggsave(snakemake@output$dispersion_plot, g, dpi = 300, height = 4, width = 6)
