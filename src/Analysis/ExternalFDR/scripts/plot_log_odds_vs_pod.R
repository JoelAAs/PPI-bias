library(ggplot2)
library(magrittr)
library(tidyverse)
library(arrow)

pod <- open_dataset(snakemake@input$pod_file) %>%
  select(pair_id, p, n_observed) %>%
  filter(n_observed > 0)

posterior <- open_dataset(snakemake@input$pair_posterior) %>%
  select(pair_id, log_odds, mean_fdr_per_test)

observed <- pod %>%
  inner_join(posterior, by = "pair_id") %>%
  mutate(log_p = log10(p)) %>%
  select(log_p, log_odds, mean_fdr_per_test) %>%
  collect()

p_log_odds_vs_pod <- ggplot(
  observed,
  aes(x = log_p, y = log_odds, color = mean_fdr_per_test)
) +
  geom_point(size = 0.3, alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(
    x = expression(log[10]("POD")),
    y = "Per-pair log odds",
    title = "Posterior log odds against probability of detection",
    color="Study FDR"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10))

ggsave(snakemake@output$plot,
       p_log_odds_vs_pod, width = 5, height = 3, dpi = 300)

## ---- never observed pairs ----------------------------------------------
# Too many to collect one row per pair, and p only varies with n_tested when
# n_observed == 0, so bin the odds and carry the multiplicity as point size.
BIN_ODDS <- 0.05

not_observed <- open_dataset(snakemake@input$pod_file) %>%
  select(pair_id, p, n_observed) %>%
  filter(n_observed == 0) %>%
  inner_join(posterior, by = "pair_id") %>%
  mutate(
    log_p = round(log10(p), 4),
    log_odds_bin = round(log_odds / BIN_ODDS) * BIN_ODDS
  ) %>%
  count(log_p, log_odds_bin, name = "count") %>%
  collect()

p_log_odds_vs_pod_not_observed <- ggplot(
  not_observed,
  aes(x = log_p, y = log_odds_bin)
) +
  geom_point(aes(size = count), colour = "darkorange", alpha = 0.6) +
  scale_size_continuous(trans = "log10") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(
    x = expression(log[10]("POD")),
    y = "Per-pair log odds",
    size = "N protein pairs",
    title = "Posterior log odds against POD, never observed pairs"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10))

ggsave(snakemake@output$plot_not_observed,
       p_log_odds_vs_pod_not_observed, width = 5, height = 3, dpi = 300)
