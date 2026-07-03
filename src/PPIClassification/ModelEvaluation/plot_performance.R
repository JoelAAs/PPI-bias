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
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("manual_figures/ROC_auc.png", g, height = 5, width = 6)

