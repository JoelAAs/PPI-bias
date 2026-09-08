library(ggplot2)
library(dplyr)
library(tidyr)

# Rank-sum AUC: how well `score` ranks the known proteins above the rest.
auc <- function(score, known) {
  r <- rank(score)
  n_pos <- sum(known)
  n_neg <- sum(!known)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  (sum(r[known]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

read_ranked <- function(path, known_ids, dataset, role) {
  df <- read.table(path, sep = "\t", header = TRUE)
  names(df)[1] <- "id"
  df %>%
    mutate(
      known = id %in% known_ids,
      dataset = dataset,
      role = role,
      facet = sprintf("%s (%s)", toupper(dataset), role)
    )
}

known_aa <- read.table(snakemake@input$known_auto_activators, sep = "\t", header = TRUE)$UniProt

crapome <- read.table(snakemake@input$crapome_proteins, sep = "\t", header = TRUE)
known_crap <- crapome$uniprot_id[
  crapome$mean_sc >= quantile(crapome$mean_sc, 1 - snakemake@params$crapome_top)
]

df <- bind_rows(
  read_ranked(snakemake@input$ranked_y2h, known_aa, "y2h", "bait"),
  read_ranked(snakemake@input$ranked_ms, known_crap, "ms", "prey")
)

g <- ggplot(df, aes(x = D, y = p)) +
  geom_point(data = ~ filter(.x, !known), alpha = 0.25, size = 0.9, color = "darkorange") +
  geom_point(data = ~ filter(.x, known), alpha = 0.9, size = 1.6, color = "firebrick") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey30") +
  facet_wrap(~facet, ncol = 2, scales = "free_x") +
  scale_y_log10() +
  labs(
    x = "Partner effect D",
    y = expression(log[10] * "(Detection rate " * italic(p) * ")"),
  ) +
  theme_bw()

ggsave(snakemake@output$partner_effect_plot, g, dpi = 300, height = 4, width = 9)

# Validation: where do the known proteins fall in each ranking?
ranked_long <- df %>%
  group_by(facet) %>%
  mutate(
    `-D` = percent_rank(-D),
    p = percent_rank(p),
    score = percent_rank(score)
  ) %>%
  ungroup() %>%
  pivot_longer(c(`-D`, p, score), names_to = "statistic", values_to = "pct_rank") %>%
  mutate(
    statistic = factor(statistic, levels = c("-D", "p", "score")),
    group = ifelse(known, "known", "other")
  )

g_val <- ggplot(ranked_long, aes(x = pct_rank, color = group)) +
  stat_ecdf(linewidth = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  facet_grid(facet ~ statistic) +
  scale_color_manual(values = c(known = "firebrick", other = "grey50")) +
  labs(
    x = "Percentile rank (1 = top of screen)",
    y = "ECDF",
    color = NULL,
    title = "Known proteins vs. background, by ranking statistic"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(snakemake@output$validation_plot, g_val, dpi = 300, height = 5, width = 9)
