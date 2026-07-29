library(ggplot2)
# library(ggbeeswarm)
library(tidyverse)
# > colnames(auc_data)
#  [1] "model"     "roc_auc"   "samples"   "dataset"
#  [5] "neg_limit" "pos_limit" "random"
#
auc_data <- read.csv(
  "/mnt/ghost/ieo7513/work_folder/classification/xgboost/permuted/all_metrics_undirectional_ESM2.csv",
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

auc_data$test_random <- sapply(
  auc_data$model, function(x) grepl("_no", x)
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

# AUC for HCNI and non-observed negative sets, plotted separately as a
# function of the positive threshold.
auc_summary <- auc_data %>%
  group_by(dataset, neg_limit, pos_limit, random) %>%
  summarise(
    n = n(),
    mean_auc = mean(roc_auc),
    se_auc = sd(roc_auc) / sqrt(n),
    ci_lo = mean_auc - qt(0.975, df = n - 1) * se_auc,
    ci_hi = mean_auc + qt(0.975, df = n - 1) * se_auc,
    .groups = "drop"
  )

g_delta <- ggplot(
  auc_summary,
  aes(
    x = pos_limit, y = mean_auc,
    color = dataset, linetype = random, group = interaction(dataset, random)
  )
) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_errorbar(
    aes(ymin = ci_lo, ymax = ci_hi),
    width = 0.15, position = position_dodge(width = 0.3)
  ) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_point(size = 2, position = position_dodge(width = 0.3)) +
  scale_color_manual(
    values = c("Combined" = "blue", "MS" = "darkgreen", "Y2H" = "darkorange")
  ) +
  scale_linetype_manual(
    values = c("FALSE" = "solid", "TRUE" = "dashed"),
    labels = c("HCNI", "Non-observed")
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
    title = "AUC over thresholds, HCNI vs. non-observed negatives",
    x = "Positive threshold",
    y = "AUC",
    color = "Dataset",
    linetype = "Negative data"
  ) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("manual_figures/ROC_auc_by_negtype.png", g_delta, height = 4, width = 6)

# Accuracy barplot: acc_I vs. acc_NI split by which negative set was used at
# test time, dodged by train_random_negative (`random`). acc_I does not
# depend on the test negative set, so it is sourced once (!test_random) to
# avoid double-counting it under both the non-observed and HRNI bars.
acc_long <- bind_rows(
  auc_data %>%
    filter(!test_random) %>%
    transmute(
      permutation, dataset, neg_limit, pos_limit, random,
      metric = "Interactions", value = acc_I
    ),
  auc_data %>%
    filter(test_random) %>%
    transmute(
      permutation, dataset, neg_limit, pos_limit, random,
      metric = "Non-interaction\n(Non-observed)", value = acc_NI
    ),
  auc_data %>%
    filter(!test_random) %>%
    transmute(
      permutation, dataset, neg_limit, pos_limit, random,
      metric = "Non-interaction\n(HRNI)", value = acc_NI
    )
)

acc_long$metric <- factor(
  acc_long$metric,
  levels = c(
    "Interactions",
    "Non-interaction\n(Non-observed)",
    "Non-interaction\n(HRNI)"
  )
)

acc_summary <- acc_long %>%
  group_by(dataset, metric, random) %>%
  summarise(
    n = n(),
    mean_acc = mean(value),
    se_acc = sd(value) / sqrt(n),
    ci_lo = mean_acc - qt(0.975, df = n - 1) * se_acc,
    ci_hi = mean_acc + qt(0.975, df = n - 1) * se_acc,
    .groups = "drop"
  )

g_acc <- ggplot(
  acc_summary,
  aes(x = metric, y = mean_acc, fill = random)
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(
    aes(ymin = ci_lo, ymax = ci_hi),
    position = position_dodge(width = 0.7), width = 0.15
  ) +
  scale_fill_manual(
    values = c("FALSE" = "darkorange", "TRUE" = "blue"),
    labels = c("HCNI", "Non-observed")
  ) +
  facet_wrap(~dataset) +
  labs(
    title = "Accuracy on interactions vs. non-interactions",
    x = NULL,
    y = "Accuracy",
    fill = "Train negative data"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = -20, hjust = 0, vjust = 1)
  )

ggsave("manual_figures/accuracy_barplot.png", g_acc, height = 4.5, width = 6)

