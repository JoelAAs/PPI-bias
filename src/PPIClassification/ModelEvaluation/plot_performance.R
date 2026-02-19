
library(ggplot2)

auc_data <- read.csv(
  "work_folder/per_gene/classification/randomforest/metrics/all_metrics.csv",
  sep="\t", header=TRUE)

auc_data$delta_pr_auc <- auc_data$pr_auc - auc_data$pr_auc_dummy
auc_data$delta_pr_auc_neg <- auc_data$pr_auc_neg - auc_data$pr_auc_dummy_neg
auc_data$delta_roc_auc <- auc_data$roc_auc - auc_data$roc_auc_dummy
auc_data$dataset <- sapply(auc_data$model, function(x) strsplit(x, "_")[[1]][1])
auc_data$partition <- sapply(auc_data$model, function(x) strsplit(x, "_")[[2]][1])

ggplot(auc_data,
       aes(x=delta_pr_auc, y=delta_pr_auc_neg, color=model)) +
  geom_point() +
  labs(title="Delta PR AUC vs Delta negative PR AUC ", x="Delta PR AUC", y="Delta NEG PR AUC") +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  facet_wrap(~ dataset) +
  theme(legend.position="bottom")

