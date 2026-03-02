library(ggplot2)
library(tidyverse)

auc_data <- read.csv(
  "work_folder/per_gene/classification/randomforest/metrics/all_metrics.csv",
  sep = "\t", header = TRUE
)

auc_data$delta_pr_auc <- auc_data$pr_auc - auc_data$pr_auc_base
auc_data$delta_pr_auc_neg <- auc_data$pr_auc_neg - auc_data$pr_auc_neg_base
auc_data$delta_roc_auc <- auc_data$roc_auc - auc_data$roc_auc_base
auc_data$dataset <- sapply(auc_data$model, function(x) strsplit(x, "_")[[1]][1])
auc_data$network_type <- sapply(
  auc_data$model, function(x) strsplit(x, "_")[[1]][2])
auc_data$selection <- sapply(
  auc_data$model, function(x) strsplit(x, "_")[[1]][3])
auc_data$partition <- sapply(
  auc_data$model, function(x) strsplit(x, "_")[[1]][4])

xy <- auc_data[auc_data$dataset == "goldensplit", c("delta_pr_auc", "delta_pr_auc_neg")]
auc_data <- auc_data[auc_data$dataset != "goldensplit", ]

g <- ggplot(
  auc_data %>% filter(network_type == "directional"),
  aes(
    x = delta_pr_auc,
    y = delta_pr_auc_neg
  )
) +
  geom_point(aes(color = dataset, shape = partition)) +
  geom_point(data = xy, color = "Black", size = 2) +
  labs(title = "Delta PR AUC vs Delta negative PR AUC (Directional)", x = "Delta PR AUC", y = "Delta NEG PR AUC") +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(selection ~ dataset) +
  theme(legend.position = "right")

ggsave("delta_directional_pr_auc.png", g, height = 4, width = 6)

g <- ggplot(
  auc_data %>% filter(network_type == "undirectional"),
  aes(
    x = delta_pr_auc,
    y = delta_pr_auc_neg
  )
) +
  geom_point(aes(color = dataset, shape = partition)) +
  geom_point(data = xy, color = "Black", size = 2) +
  labs(title = "Delta PR AUC vs Delta negative PR AUC (Undirected)", x = "Delta PR AUC", y = "Delta NEG PR AUC") +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(selection ~ dataset) +
  theme(legend.position = "right")


ggsave("delta_undirectional_pr_auc.png", g, height = 4, width = 6)
