library(tidyverse)
library(arrow)

id_cols <- paste0(snakemake@params$id_pattern, c("_bait", "_prey"))

pod <- open_dataset(snakemake@input$pod_file) %>%
  select(all_of(id_cols), p, lower_bound_pod)

observed <- open_dataset(snakemake@input$pair_posterior) %>%
  filter(n_observed > 0) %>%
  select(all_of(id_cols), posterior, log_odds, n_tested) %>%
  inner_join(pod, by = id_cols) %>%
  collect() %>%
  arrange(n_tested) # most-tested pairs drawn on top

p_fdr_vs_pod <- ggplot(observed, aes(x = lower_bound_pod, y = posterior)) +
  geom_point(size = 0.3, alpha = 0.2) +
  scale_x_log10() +
  labs(
    x = "POD lower bound",
    y = "Interaction probability",
    title = "Per-pair posterior against POD lower bound, observed pairs"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10))

ggsave(
  snakemake@output$fdr_vs_pod,
  p_fdr_vs_pod,
  width = 5,
  height = 3,
  dpi = 300
)

## ---- log odds, POD and n_tested --------------------------------------------
theme_pod <- theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 10)
  )
odds_zero <- geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40")
median_line <- stat_summary(
  fun = median,
  geom = "line",
  colour = "#d95f02",
  linewidth = 0.7
)

p_log_odds_vs_pod <- ggplot(
  observed,
  aes(x = p, y = log_odds, colour = n_tested)
) +
  geom_point(size = 0.3, alpha = 0.3) +
  odds_zero +
  scale_x_log10() +
  scale_colour_viridis_c(trans = "log10") +
  labs(
    x = "POD",
    y = "Posterior log odds",
    colour = "N tested",
    title = "Posterior log odds against POD, observed pairs"
  ) +
  theme_pod

p_log_odds_vs_n_tested <- ggplot(observed, aes(x = n_tested, y = log_odds)) +
  geom_point(size = 0.3, alpha = 0.1) +
  odds_zero +
  median_line +
  scale_x_log10() +
  labs(
    x = "N tested",
    y = "Posterior log odds",
    title = "Posterior log odds against times tested, orange: median"
  ) +
  theme_pod

p_pod_vs_n_tested <- ggplot(observed, aes(x = n_tested, y = p)) +
  geom_point(size = 0.3, alpha = 0.1) +
  median_line +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "N tested",
    y = "POD",
    title = "POD against times tested, orange: median"
  ) +
  theme_pod

ggsave(
  snakemake@output$log_odds_vs_pod,
  p_log_odds_vs_pod,
  width = 5,
  height = 3.5,
  dpi = 300
)
ggsave(
  snakemake@output$log_odds_vs_n_tested,
  p_log_odds_vs_n_tested,
  width = 5,
  height = 3.5,
  dpi = 300
)
ggsave(
  snakemake@output$pod_vs_n_tested,
  p_pod_vs_n_tested,
  width = 5,
  height = 3.5,
  dpi = 300
)

## ---- theta_s distribution --------------------------------------------------
study <- read_tsv(snakemake@input$study_fdr, show_col_types = FALSE) %>%
  mutate(
    evidence = if_else(n_informative > 0, "External evidence", "Prior only")
  )
pi_ref <- read_tsv(snakemake@input$global_fit, show_col_types = FALSE)$pi

p_theta_s <- ggplot(study, aes(x = theta_s, fill = evidence)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = pi_ref, linetype = "dashed", colour = "grey40") +
  scale_x_log10() +
  labs(
    x = expression(theta[s]),
    y = "N studies",
    fill = NULL,
    title = "Study precision, dashed line at reference pi (floor for w_detected)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10))

ggsave(snakemake@output$theta_s, p_theta_s, width = 5, height = 3)
