library(ggplot2)
library(dplyr)

read_detectability <- function(path) {
  df <- read.table(path, sep = "\t", header = TRUE)
  names(df)[1] <- "id"
  df
}

bait <- read_detectability(snakemake@input$bait_detectability)
prey <- read_detectability(snakemake@input$prey_detectability)

bait_only <- setdiff(bait$id, prey$id)
prey_only <- setdiff(prey$id, bait$id)
cat(sprintf(
  "Skipping %d unique protein(s) that only occur as bait or prey (%d bait-only, %d prey-only)\n",
  length(union(bait_only, prey_only)), length(bait_only), length(prey_only)
))

df <- inner_join(
  bait[, c("id", "log_odds_deviation", "log_odds_sd")],
  prey[, c("id", "log_odds_deviation", "log_odds_sd")],
  by = "id",
  suffix = c("_bait", "_prey")
) %>%
  mutate(max_sd = pmax(log_odds_sd_bait, log_odds_sd_prey))

localisation <- read.table(snakemake@input$protein_localisation, sep = "\t", header = TRUE)

no_localisation <- setdiff(df$id, localisation$id)
cat(sprintf(
  "Skipping %d unique protein(s) without localisation data\n",
  length(no_localisation)
))

df_localisation <- inner_join(df, localisation, by = "id")

top_localisations <- df_localisation %>%
  distinct(id, localisation) %>%
  count(localisation, sort = TRUE) %>%
  slice_head(n = 3) %>%
  pull(localisation)

df_localisation <- df_localisation %>%
  filter(localisation %in% top_localisations) %>%
  mutate(localisation = factor(localisation, levels = top_localisations))

sd_breaks <- quantile(df_localisation$max_sd, probs = c(0, 0.1, 0.3, 1))
percentile_names <- c("0-10%", "10-30%", "30-100%")

df_localisation <- df_localisation %>%
  mutate(sd_bin = cut(
    max_sd, breaks = sd_breaks,
    include.lowest = TRUE, labels = percentile_names
  ))

stat_labels <- df_localisation %>%
  group_by(localisation, sd_bin) %>%
  summarise(
    rho = cor(log_odds_deviation_bait, log_odds_deviation_prey, method = "spearman"),
    r = cor(log_odds_deviation_bait, log_odds_deviation_prey, method = "pearson"),
    slope = coef(lm(log_odds_deviation_prey ~ log_odds_deviation_bait))[2],
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf(
    "Spearman rho = %.3f\nPearson r = %.3f, slope = %.3f (n = %d)",
    rho, r, slope, n
  ))

g <- ggplot(df_localisation, aes(x = log_odds_deviation_bait, y = log_odds_deviation_prey)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "grey40", linetype = "dashed") +
  geom_text(
    data = stat_labels,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.05, vjust = 1.2, inherit.aes = FALSE, size = 2.5
  ) +
  facet_grid(localisation ~ sd_bin) +
  labs(x = expression(b[i]), y = expression(p[i])) +
  theme_bw()

ggsave(snakemake@output$plot, g, dpi = 300, height = 10, width = 10)
