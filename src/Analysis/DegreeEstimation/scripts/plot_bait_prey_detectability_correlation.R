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
  bait[, c("id", "detectability")],
  prey[, c("id", "detectability")],
  by = "id",
  suffix = c("_bait", "_prey")
) %>%
  mutate(
    log_detectability_bait = log(detectability_bait),
    log_detectability_prey = log(detectability_prey)
  )

stat_label <- df %>%
  summarise(
    rho = cor(log_detectability_bait, log_detectability_prey, method = "spearman"),
    r = cor(log_detectability_bait, log_detectability_prey, method = "pearson"),
    slope = coef(lm(log_detectability_prey ~ log_detectability_bait))[2],
    n = n()
  ) %>%
  mutate(label = sprintf(
    "Spearman rho = %.3f\nPearson r = %.3f, slope = %.3f (n = %d)",
    rho, r, slope, n
  ))

g <- ggplot(df, aes(x = log_detectability_bait, y = log_detectability_prey)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "grey40", linetype = "dashed") +
  geom_text(
    data = stat_label,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.05, vjust = 1.2, inherit.aes = FALSE, size = 3
  ) +
  labs(x = "log(bait detectability)", y = "log(prey detectability)") +
  theme_bw()

ggsave(snakemake@output$plot, g, dpi = 300, height = 5, width = 5)
