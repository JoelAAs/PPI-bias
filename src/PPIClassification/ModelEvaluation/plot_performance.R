library(ggplot2)
#library(ggbeeswarm)
library(tidyverse)
# > colnames(auc_data)
#  [1] "model"           "pr_auc"          "pr_auc_base"     "pr_auc_neg"     
#  [5] "pr_auc_neg_base" "roc_auc"         "roc_auc_base"    "ce_obs"         
#  [9] "ce_baseline"     "samples"         "dataset"         "neg_limit"      
# [13] "pos_limit"       "random"         
# > 
auc_data <- read.csv(
  "work_folder/per_gene/classification/xgboost/metrics/all_metrics.csv",
  sep = "\t", header = TRUE
)

auc_data$dataset <- sapply(auc_data$model, function(x) strsplit(x, "_")[[1]][1])
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
  levels = c(1, 2)   # desired order
)

auc_data$pos_limit <- factor(
  auc_data$pos_limit,
  levels = c("all", "0.02", "0.15")  # desired order
)



g <- ggplot(
  auc_data,
  aes(
    y = pr_auc,
    x = pr_auc_neg
  )
) +
  geom_point(aes(color = random, shape = dataset)) +
  labs(
    title = "PR AUC for interaction and non-interaction classification",
    x = "Non-interaction PR AUC",
    y = "Interaction PR AUC",
    color = "Negative data",
    shape = "Data type"
  ) +
  scale_color_manual(
    values = c("darkorange", "blue"),
    labels = c("HCNI", "Non-observed")
  ) +
  scale_shape_manual(
    values = c(15,16,17),
    labels = c("Combined", "MS", "Y2H")
  ) +
  theme_bw() +
  facet_grid(
    neg_limit ~pos_limit,
    labeller = labeller(
      neg_limit = c(
        "1" = "Negative tests >= 1",
        "2" = "Negative tests >= 2"
      ),
      pos_limit = c(
        "all" = "Any interactions",
        "0.02" = "P- >= 0.02",
        "0.15" = "P- >= 0.15"
      )
    )
  ) +  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_vline(xintercept = 0.95, linetype = "dashed") +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("manual_figures/PR_auc.png", g, height = 4, width = 6)



g <- ggplot(
  auc_data,
  aes(
    x = dataset,
    y = ce_obs
  )
) +
  geom_point(aes(color = random)) +
  labs(
    title = "Cross entropy between negative data selection strategies",
    y = "Cross-entropy",
    x = "Dataset",
    color = "Negative data"
  ) +
  scale_color_manual(
    values = c("darkorange", "blue"),
    labels = c("HCNI", "Non-observed")
  ) +
  theme_bw() +
  scale_x_discrete(
  labels = c(
    "flat" = "Combined",
    "ms" = "MS",
    "y2h" = "Y2H"
  )
) +
  facet_grid(
    neg_limit ~pos_limit,
    labeller = labeller(
      neg_limit = c(
        "1" = "Negative tests >= 1",
        "2" = "Negative tests >= 2"
      ),
      pos_limit = c(
        "all" = "Any interactions",
        "0.02" = "P.025 >= 0.02",
        "0.15" = "P.025 >= 0.15"
      )
    )
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("manual_figures/CE.png", g, height = 4, width = 6)



degree_balance <- read.table("work_folder/per_gene/subsets/train/equal_edge/balance/all_metrics.csv", sep="\t", header=T)

degree_balance$total_delta <- degree_balance$bait_degree_delta + degree_balance$prey_degree_delta
degree_balance$neg_type <- ifelse(degree_balance$random, "Non-observed", "HCNI")
degree_balance$dataset <- factor(
  degree_balance$dataset,
  labels = c("flat" = "Combined", "ms" = "MS", "y2h" = "Y2H")
)

g_degree <- ggplot(
  degree_balance,
  aes(x = dataset, y = total_delta/num_edges, fill = neg_type)
) +
  geom_boxplot(outlier.size = 0.8) +
  labs(
    title = "Train set degree balance by dataset and \nnegative sampling strategy",
    x = "Dataset",
    y = "Degree delta per edge (bait + prey)/(|E+|+ |E-|)",
    fill = "Negative data"
  ) +
  scale_fill_manual(values = c("HCNI" = "darkorange", "Non-observed" = "blue")) +
  theme_bw()

ggsave("manual_figures/degree_balance.png", g_degree, height = 5, width = 6)

g_edges <- ggplot(
  degree_balance[degree_balance$random != "True", ],
  aes(x = dataset, y = num_edges)
) +
  geom_boxplot(outlier.size = 0.8, fill = "darkorange") +
  labs(
    title = "Number of edges per train set",
    x = "Dataset",
    y = "Number of edges (|E+| + |E-|)"
  ) +
  theme_bw()

ggsave("manual_figures/num_edges.png", g_edges, height = 5, width = 6)

