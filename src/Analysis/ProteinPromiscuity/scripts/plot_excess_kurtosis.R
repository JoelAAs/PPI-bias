library(ggplot2)
library(dplyr)

# A log-normal detectability is normal on log10, i.e. zero excess kurtosis. A heavy tail
# (a subpopulation of indiscriminately detected proteins) shows up as positive excess.

read_adjusted <- function(path, dataset) {
  df <- read.table(path, sep = "\t", header = TRUE)
  data.frame(
    dataset = dataset,
    log_detectability = log10(df$adjusted_detectability),
    min_n_pairs = df$min_n_pairs
  )
}

# D'Agostino-style test: bias-corrected excess kurtosis over its standard error under
# normality. At these n the normal approximation is close enough that the Anscombe-Glynn
# refinement changes nothing material.
kurtosis_test <- function(x) {
  n <- length(x)
  m2 <- mean((x - mean(x))^2)
  m4 <- mean((x - mean(x))^4)
  g2 <- m4 / m2^2 - 3
  G2 <- ((n + 1) * g2 + 6) * (n - 1) / ((n - 2) * (n - 3))
  se <- sqrt(24 * n * (n - 1)^2 / ((n - 3) * (n - 2) * (n + 3) * (n + 5)))
  data.frame(excess_kurtosis = G2, z = G2 / se, p_value = 2 * pnorm(-abs(G2 / se)), n = n)
}

skewness <- function(x) {
  mean((x - mean(x))^3) / mean((x - mean(x))^2)^1.5
}

min_pairs <- snakemake@params$min_pairs

filtered_label <- sprintf("min_n_pairs >= %d", min_pairs)

all_proteins <- bind_rows(
  read_adjusted(snakemake@input$ms_adjusted_detectability, "ms"),
  read_adjusted(snakemake@input$y2h_adjusted_detectability, "y2h")
)
# the tail is quantified on the well-tested proteins only; the unfiltered set is drawn
# alongside to show what the min_n_pairs cut removes, but is not fitted or tested
df <- filter(all_proteins, min_n_pairs >= min_pairs)

for (ds in unique(all_proteins$dataset)) {
  cat(sprintf(
    "%s: %d of %d protein(s) with %s\n", ds,
    sum(df$dataset == ds), sum(all_proteins$dataset == ds), filtered_label
  ))
}

densities <- bind_rows(
  mutate(all_proteins, subset = "all"),
  mutate(df, subset = filtered_label)
) %>%
  mutate(subset = factor(subset, levels = c(filtered_label, "all")))

stats <- bind_rows(lapply(split(df, df$dataset), function(d) {
  cbind(dataset = d$dataset[1], kurtosis_test(d$log_detectability),
        skewness = skewness(d$log_detectability))
}))
print(stats, row.names = FALSE)

stat_labels <- stats %>%
  left_join(count(all_proteins, dataset, name = "n_all"), by = "dataset") %>%
  mutate(label = sprintf(
    "excess kurtosis = %+.3f\n n = %d ",
    excess_kurtosis, n
  ))

# log-normal reference: the normal fitted to each dataset's own log10 detectability
normal_ref <- bind_rows(lapply(split(df, df$dataset), function(d) {
  grid <- seq(min(d$log_detectability), max(d$log_detectability), length.out = 400)
  data.frame(
    dataset = d$dataset[1],
    log_detectability = grid,
    density = dnorm(grid, mean(d$log_detectability), sd(d$log_detectability))
  )
}))

g <- ggplot(densities, aes(x = log_detectability)) +
  geom_density(aes(color = dataset, linetype = subset), linewidth = 0.8,
               key_glyph = "path") +
  geom_line(
    data = normal_ref, aes(y = density),
    linetype = "dashed", color = "grey40"
  ) +
  geom_label(
    data = stat_labels,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.05, vjust = 1.2, inherit.aes = FALSE, size = 3,
    label.size = 0.3, fill = "white"
  ) +
  facet_wrap(~dataset, ncol = 2, labeller = as_labeller(toupper), scales = "free_y") +
  scale_color_manual(values = c(ms = "forestgreen", y2h = "darkorange"), guide = "none") +
  scale_linetype_manual(values = setNames(c("solid", "dotted"), c(filtered_label, "all"))) +
  labs(
    x = expression(log[10] ~ p[avg]),
    y = "Density",
    linetype = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(snakemake@output$plot_distribtuions, g, dpi = 300, height = 4, width = 8)
