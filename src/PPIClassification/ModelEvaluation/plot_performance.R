library(ggplot2)
library(tidyverse)

auc_data <- read.csv(
  "work_folder/per_gene/classification/randomforest/metrics/all_metrics.csv",
  sep = "\t", header = TRUE
)

auc_data$delta_pr_auc <- auc_data$pr_auc - auc_data$pr_auc_base
auc_data$delta_pr_auc_neg <- auc_data$pr_auc_neg - auc_data$pr_auc_neg_base
auc_data$delta_roc_auc <- auc_data$roc_auc - auc_data$roc_auc_base
auc_data$delta_ce <- auc_data$ce_obs - auc_data$ce_baseline
auc_data$dataset <- sapply(auc_data$model, function(x) strsplit(x, "_")[[1]][1])
auc_data$network_type <- sapply(
  auc_data$model, function(x) strsplit(x, "_")[[1]][2]
)
auc_data$selection <- sapply(
  auc_data$model, function(x) strsplit(x, "_")[[1]][3]
)
auc_data$partition <- sapply(
  auc_data$model, function(x) strsplit(strsplit(x, "_")[[1]][4], "-")[[1]][1]
)
auc_data$random <- sapply(
  auc_data$model, function(x) grepl("random", x)
)


mat <- auc_data %>% 
  filter(random == TRUE & partition == "sequencesimilarity" & network_type == "directional") %>%
  {.[, c("dataset", "samples", "selection")]} %>%
  pivot_wider(
    names_from = dataset,
    values_from = samples
  )
  column_to_rownames("row_id") %>%
  as.matrix()


xy <- auc_data[auc_data$dataset == "goldensplit", c("delta_pr_auc", "delta_pr_auc_neg")]
auc_data <- auc_data[auc_data$dataset != "goldensplit", ]

g <- ggplot(
  auc_data, 
  aes(
    y = -delta_ce,
    x = interaction(partition, network_type)
  )
) +
  geom_point(aes(color = random, shape = network_type)) +
  labs(
    title = expression(paste(Delta, " ", Delta, "  Log-loss between HCNI and non-observed")),
    x = "Gene partition",
    y = expression(paste(Delta, " ", Delta, "  Log-loss")),
    color = "Negative data",
    shape = "Network type"
  ) +
  scale_color_manual(
    values = c("darkorange", "blue"),
    labels = c("HCNI", "Non-observed")
  ) +
  theme_bw() +
  facet_grid(selection ~ dataset) +
  scale_x_discrete(labels = function(x) sub("\\..*", "", x)) + 
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("delta_directional_ce.png", g, height = 4, width = 6)

g <- ggplot(
  auc_data,
  aes(
    y = roc_auc,
    x = interaction(partition, network_type)
 )) +
  geom_point(aes(color = random, shape = network_type)) +
  labs(
    title = "ROC AUC between HCNI and non-observed",
    x = "Gene partition",
    y = "ROC AUC",
    color = "Negative data",
    shape = "Network type"
  ) +
  scale_color_manual(
    values = c("darkorange", "blue"),
    labels = c("HCNI", "Non-observed")
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_bw() +
  facet_grid(selection ~ dataset) +
  scale_x_discrete(labels = function(x) sub("\\..*", "", x)) + 
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("model_roc.png", g, height = 4, width = 6)


df_balance_undir <- read.csv(
  "work_folder/per_gene/subsets/degree_balance/all_undirectional.csv",
  sep = "\t", header = TRUE
)

df_balance_undir$model <- df_balance_undir$dataset
df_balance_undir$dataset <- sapply(df_balance_undir$model, function(x) strsplit(x, "_")[[1]][1])
df_balance_undir$negative_lim <- sapply(
  df_balance_undir$model, function(x) strsplit(x, "_")[[1]][4]
)
df_balance_undir$positive_lim <- sapply(
  df_balance_undir$model, function(x) strsplit(x, "_")[[1]][6]
)
df_balance_undir$partition <- sapply(
  df_balance_undir$model, function(x) strsplit(strsplit(x, "_")[[1]][7], "-")[[1]][1]
)
df_balance_undir$random <- sapply(
  df_balance_undir$model, function(x) grepl("random", x)
)
df_balance_undir <- df_balance_undir %>%
  mutate(selection = case_when(
    negative_lim == 1 & positive_lim == 0.02 ~ "loose",
    negative_lim == 2 & positive_lim == 0.15 ~ "medium",
    negative_lim == 3 & positive_lim == 0.29 ~ "strict",
    TRUE ~ NA_character_
  ))

g <- ggplot(
  df_balance_undir,
  aes(
    y = p_undir,
    x = sp_undir
  )
) +
  geom_point(aes(color = random, shape = set_type)) +
  labs(
    title = "Degree balance score between HCNI and non-observed\nUndirected",
    x = "Spearman correlation",
    y = "Pearson correlation",
    color = "Negative data"
  ) +
  scale_color_manual(
    values = c("darkorange", "blue"),
    labels = c("HCNI", "Non-observed")
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  facet_grid(selection ~ dataset) +
  #scale_x_discrete(labels = function(x) sub("\\..*", "", x)) + 
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("degree_balance_undir.png",g,  height = 4, width = 6)



df_balance_dir <- read.csv(
  "work_folder/per_gene/subsets/degree_balance/all_directional.csv",
  sep = "\t", header = TRUE
)

df_balance_dir$model <- df_balance_dir$dataset
df_balance_dir$dataset <- sapply(df_balance_dir$model, function(x) strsplit(x, "_")[[1]][1])
df_balance_dir$negative_lim <- sapply(
  df_balance_dir$model, function(x) strsplit(x, "_")[[1]][4]
)
df_balance_dir$positive_lim <- sapply(
  df_balance_dir$model, function(x) strsplit(x, "_")[[1]][6]
)
df_balance_dir$partition <- sapply(
  df_balance_dir$model, function(x) strsplit(strsplit(x, "_")[[1]][7], "-")[[1]][1]
)
df_balance_dir$random <- sapply(
  df_balance_dir$model, function(x) grepl("random", x)
)
df_balance_dir <- df_balance_dir %>%
  mutate(selection = case_when(
    negative_lim == 1 & positive_lim == 0.02 ~ "loose",
    negative_lim == 2 & positive_lim == 0.15 ~ "medium",
    negative_lim == 3 & positive_lim == 0.29 ~ "strict",
    TRUE ~ NA_character_
  ))

g <- ggplot(
  df_balance_dir,
  aes(
    y = sp_bait,
    x = sp_prey
  )
) +
  geom_point(aes(color = random, shape = set_type)) +
  labs(
    title = "Degree balance score between HCNI and non-observed\nDirected",
    x = "Spearman correlation (prey)",
    y = "Spearman correlation (bait)",
    color = "Negative data"
  ) +
  scale_color_manual(
    values = c("darkorange", "blue"),
    labels = c("HCNI", "Non-observed")
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  facet_grid(selection ~ dataset) +
  #scale_x_discrete(labels = function(x) sub("\\..*", "", x)) + 
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("degree_balance_dir.png",g,  height = 4, width = 6)
