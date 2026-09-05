library(ggplot2)
library(dplyr)

read_keyed <- function(path) {
  df <- read.table(path, sep = "\t", header = TRUE)
  names(df)[1] <- "id"
  df
}

df <- inner_join(
  read_keyed(snakemake@input$ms_sticky_proteins),
  read_keyed(snakemake@input$ms_prey_detectability),
  by = "id"
)

stat_label <- sprintf(
  "Spearman rho = %.3f (n = %d)",
  cor(df$detectability, df$dispersion_ratio, method = "spearman"),
  nrow(df)
)

# a prey detected at the same rate by every bait has obs_var == 0; floor those
# onto the log axis rather than dropping them, they are the sticky ones
floor_ratio <- min(df$dispersion_ratio[df$dispersion_ratio > 0]) / 2
df$dispersion_ratio_floored <- pmax(df$dispersion_ratio, floor_ratio)

g <- ggplot(df, aes(x = detectability, y = dispersion_ratio_floored)) +
  geom_point(alpha = 0.4, size = 1, color = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  annotate(
    "label", x = -Inf, y = Inf, label = stat_label,
    hjust = -0.05, vjust = 1.2, size = 3, linewidth = 0.3, fill = "white"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "log(MS prey detectability)",
    y = "log(Dispersion ratio: observed / expected variance)"
  ) +
  theme_bw()

ggsave(snakemake@output$sticky_detectable_plot, g, dpi = 300, height = 4, width = 6)
