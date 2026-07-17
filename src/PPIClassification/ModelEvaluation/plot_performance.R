library(ggplot2)
# library(ggbeeswarm)
library(tidyverse)
# > colnames(auc_data)
#  [1] "model"     "roc_auc"   "samples"   "dataset"
#  [5] "neg_limit" "pos_limit" "random"
# >
auc_data <- read.csv(
  "work_folder/classification/xgboost/permuted/all_metrics_undirectional_ESM2.csv",
  sep = "\t", header = TRUE
)

auc_data$dataset <- sapply(auc_data$model, function(x) strsplit(x, "_")[[1]][1])
auc_data$dataset <- sapply(auc_data$dataset, function(x) {
  case_when(
    x == "flat" ~ "Combined",
    x == "ms" ~ "MS",
    x == "y2h" ~ "Y2H",
    TRUE ~ as.character(x)
  )
})
auc_data$neg_limit <- sapply(
  auc_data$model, function(x) strsplit(x, "_")[[1]][4]
)
auc_data$pos_limit <- sapply(
  auc_data$model, function(x) strsplit(strsplit(x, "_")[[1]][6], "-")[[1]][1]
)
auc_data$random <- sapply(
  auc_data$model, function(x) grepl("random", x)
)
auc_data$neg_limit <- factor(
  auc_data$neg_limit,
  levels = c(1, 2) # desired order
)

auc_data$pos_limit <- factor(
  auc_data$pos_limit,
  levels = c("all", "0.02", "0.15") # desired order
)

# plotmath expressions shared by both plots' pos_limit/neg_limit labels
pos_limit_labels <- c(
  "all" = "Any~interactions",
  "0.02" = "Q[2.5] > 0.02",
  "0.15" = "Q[2.5] > 0.15"
)

neg_limit_labels <- c(
  "1" = "Negative~tests >= 1",
  "2" = "Negative~tests >= 2"
)


g <- ggplot(
  auc_data,
  aes(
    x = dataset,
    y = roc_auc
  )
) +
  geom_boxplot(aes(color = random)) +
  labs(
    title = "ROC AUC per datasets and threshold configurations",
    x = "Detection dataset",
    y = "ROC AUC",
    color = "Negative data",
    shape = "Data type"
  ) +
  scale_color_manual(
    values = c("darkorange", "blue"),
    labels = c("HCNI", "Non-observed")
  ) +
  theme_bw() +
  facet_grid(
    neg_limit ~ pos_limit,
    labeller = labeller(
      neg_limit = as_labeller(neg_limit_labels, label_parsed),
      pos_limit = as_labeller(pos_limit_labels, label_parsed)
    )
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("manual_figures/ROC_auc.png", g, height = 5, width = 6)

# ΔAUC = AUC_HCNI - AUC_nonobs per replicate, to show the HCNI advantage
# directly as a function of the positive threshold.
hcni_df <- auc_data %>%
  filter(!random) %>%
  select(permutation, dataset, neg_limit, pos_limit, roc_auc_hcni = roc_auc)

nonobs_df <- auc_data %>%
  filter(random) %>%
  select(permutation, dataset, neg_limit, pos_limit, roc_auc_nonobs = roc_auc)

delta_df <- inner_join(
  hcni_df, nonobs_df,
  by = c("permutation", "dataset", "neg_limit", "pos_limit")
) %>%
  mutate(delta_auc = roc_auc_hcni - roc_auc_nonobs)

delta_summary <- delta_df %>%
  group_by(dataset, neg_limit, pos_limit) %>%
  summarise(
    n = n(),
    mean_delta = mean(delta_auc),
    se_delta = sd(delta_auc) / sqrt(n),
    ci_lo = mean_delta - qt(0.975, df = n - 1) * se_delta,
    ci_hi = mean_delta + qt(0.975, df = n - 1) * se_delta,
    .groups = "drop"
  )

g_delta <- ggplot(
  delta_summary,
  aes(
    x = pos_limit, y = mean_delta,
    color = dataset, group = dataset
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbar(
    aes(ymin = ci_lo, ymax = ci_hi),
    width = 0.15, position = position_dodge(width = 0.3)
  ) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_point(size = 2, position = position_dodge(width = 0.3)) +
  scale_color_manual(
    values = c("Combined" = "blue", "MS" = "darkgreen", "Y2H" = "darkorange")
  ) +
  scale_x_discrete(
    labels = parse(text = pos_limit_labels)
  ) +
  facet_wrap(
    ~ neg_limit,
    labeller = labeller(
      neg_limit = as_labeller(neg_limit_labels, label_parsed)
    )
  ) +
  labs(
    title = "HCNI advantage over non-observed negatives",
    x = "Positive threshold",
    y = expression(Delta * "AUC (HCNI - Non-observed)"),
    color = "Dataset"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("manual_figures/ROC_auc_delta.png", g_delta, height = 4, width = 6)

