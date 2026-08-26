library(ggplot2)
library(tidyverse)

# columns: dataset, network_type, pos_limit, neg_limit, permutation, label, jaccard, n_pairs
jaccard_data <- read.csv(
  "/mnt/ghost/ieo7513/work_folder/classification/xgboost/full_test_predictions/jaccard/all_jaccard_ESM2_undirectional_similarity.tsv",
  sep = "\t", header = TRUE
)

jaccard_data$pos_limit <- factor(jaccard_data$pos_limit)
jaccard_data$neg_limit <- factor(jaccard_data$neg_limit)

# mean/CI across permutations, for each dataset x label x threshold combo
jaccard_summary <- jaccard_data %>%
  group_by(dataset, label, pos_limit, neg_limit) %>%
  summarise(
    n = n(),
    mean_jaccard = mean(jaccard, na.rm = TRUE),
    se_jaccard = sd(jaccard, na.rm = TRUE) / sqrt(n),
    ci_lo = mean_jaccard - qt(0.975, df = n - 1) * se_jaccard,
    ci_hi = mean_jaccard + qt(0.975, df = n - 1) * se_jaccard,
    .groups = "drop"
  )

g_pos <- ggplot(
  jaccard_summary,
  aes(x = pos_limit, y = mean_jaccard, color = dataset, group = dataset)
) +
  geom_errorbar(
    aes(ymin = ci_lo, ymax = ci_hi),
    width = 0.15, position = position_dodge(width = 0.3)
  ) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_point(size = 2, position = position_dodge(width = 0.3)) +
  facet_grid(label ~ neg_limit) +
  labs(
    title = "HRNI vs. random-negative model agreement by positive threshold",
    x = "Positive threshold",
    y = "Jaccard index (correct-prediction overlap)",
    color = "Dataset"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("manual_figures/jaccard_vs_pos_limit.png", g_pos, height = 7, width = 8)

g_neg <- ggplot(
  jaccard_summary,
  aes(x = neg_limit, y = mean_jaccard, color = dataset, group = dataset)
) +
  geom_errorbar(
    aes(ymin = ci_lo, ymax = ci_hi),
    width = 0.15, position = position_dodge(width = 0.3)
  ) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_point(size = 2, position = position_dodge(width = 0.3)) +
  facet_grid(label ~ pos_limit) +
  labs(
    title = "HRNI vs. random-negative model agreement by negative threshold",
    x = "Negative threshold",
    y = "Jaccard index (correct-prediction overlap)",
    color = "Dataset"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("manual_figures/jaccard_vs_neg_limit.png", g_neg, height = 7, width = 8)
