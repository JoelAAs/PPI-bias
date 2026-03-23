library(ggplot2)
library(ggbeeswarm)
library(tidyverse)

auc_data <- read.csv(
  "work_folder/per_gene/classification/randomforest/metrics/all_metrics.csv",
  sep = "\t", header = TRUE
)

auc_data$delta_pr_auc <- auc_data$pr_auc - auc_data$pr_auc_base
auc_data$delta_pr_auc_neg <- auc_data$pr_auc_neg - auc_data$pr_auc_neg_base
auc_data$delta_roc_auc <- auc_data$roc_auc - auc_data$roc_auc_base
auc_data$delta_ce <- auc_data$ce_baseline - auc_data$ce_obs
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
  {
    .[, c("dataset", "samples", "selection")]
  } %>%
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
    y = delta_ce,
    x = interaction(partition, network_type)
  )
) +
  geom_point(aes(color = random, shape = network_type)) +
  labs(
    title = expression(paste(H[baseline], " - ", H[obs], " between HCNI and non-observed")),
    x = "Gene partition",
    y = expression(paste(H[baseline], " - ", H[obs])),
    color = "Negative data",
    shape = "Network type"
  ) +
  scale_color_manual(
    values = c("darkorange", "blue"),
    labels = c("HCNI", "Non-observed")
  ) +
  theme_bw() +
  facet_grid(selection ~ dataset, labeller = labeller(
    dataset = c(
      flat = "Combined",
      ms = "MS",
      y2h = "Y2H"
    )
  )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(labels = function(x) {
    case_when(
      str_remove(x, "\\..*") == "maxpos" ~ "MaxPos",
      str_remove(x, "\\..*") == "sequencesimilarity" ~ "SeqSim",
      TRUE ~ x
    )
  }) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("delta_ce.png", g, height = 4, width = 6)

g <- ggplot(
  auc_data,
  aes(
    x = delta_pr_auc / (1 - pr_auc_base),
    y = delta_pr_auc_neg / (1 - pr_auc_neg_base)
  )
) +
  geom_point(aes(color = interaction(random, partition), shape = network_type)) +
  labs(
    title = "PR AUC between HCNI and non-observed",
    x = "(PR AUC - p)/(1-p) (interaction)",
    y = "(PR AUC - p)/(1-p) (non-interaction)",
    color = "Negative data",
    shape = "Network type"
  ) +
  scale_color_manual(
    values = c("darkorange", "blue", "#BA3F1D", "#80A1C1"),
    labels = c(
      "HCNI:MaxPos", "Non-observed:MaxPos",
      "HCNI:SeqSim", "Non-observed:SeqSim"
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  facet_grid(selection ~ dataset, labeller = labeller(
    dataset = c(
      flat = "Combined",
      ms = "MS",
      y2h = "Y2H"
    )
  )) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("model_pr_auc.png", g, height = 4, width = 6)


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
    y = sp_undir,
    x = partition
  )
) +
  geom_beeswarm(aes(color = random, shape = set_type), cex = 3, priority = "density", size = 1) +
  labs(
    title = "Correlation of Positive and negative degree:\nInteraction network",
    x = "Gene partition",
    y = "Spearman correlation",
    color = "Negative data",
    shape = "Set"
  ) +
  scale_color_manual(
    values = c("darkorange", "blue"),
    labels = c("HCNI", "Non-observed")
  ) +
  ylim(0, 1) +
  theme_bw() +
  facet_grid(selection ~ dataset, labeller = labeller(
    dataset = c(
      flat = "Combined",
      ms = "MS",
      y2h = "Y2H"
    )
  )) +
  scale_x_discrete(labels = function(x) {
    case_when(
      str_remove(x, "\\..*") == "maxpos" ~ "MaxPos",
      str_remove(x, "\\..*") == "sequencesimilarity" ~ "SeqSim",
      TRUE ~ x
    )
  }) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("degree_balance_undir.png", g, height = 4, width = 6)


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
  geom_point(aes(color = interaction(random, partition), shape = set_type)) +
  labs(
    title = "Correlation of Positive and negative degree:\nBait-prey network",
    x = "Spearman correlation (prey)",
    y = "Spearman correlation (bait)",
    color = "Negative data",
    shape = "Set"
  ) +
  scale_color_manual(
    values = c("darkorange", "blue", "#BA3F1D", "#80A1C1"),
    labels = c(
      "HCNI:MaxPos", "Non-observed:MaxPos",
      "HCNI:SeqSim", "Non-observed:SeqSim"
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  facet_grid(selection ~ dataset, labeller = labeller(
    dataset = c(
      flat = "Combined",
      ms = "MS",
      y2h = "Y2H"
    )
  )) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)
  )

ggsave("degree_balance_dir.png", g, height = 4, width = 6)

df_edge_count <- read.table(
  "edge_count.csv",
  header = F, col.names = c("model", "n_edges")
)

df_edge_count$dataset <- sapply(
  df_edge_count$model, function(x) strsplit(x, "_")[[1]][1]
)
df_edge_count$positive_lim <- sapply(
  df_edge_count$model, function(x) strsplit(x, "_")[[1]][4]
)

df_edge_count$network_type <- sapply(
  df_edge_count$model, function(x) strsplit(x, "_")[[1]][2]
)

df_edge_count_dir <- df_edge_count %>%
  filter(network_type == "directional")

df_joined <- df_balance_dir %>%
  filter(random == FALSE) %>%
  left_join(df_edge_count_dir, by = c("dataset", "positive_lim"))


df_joined[, "percent_retained"] <- 100 * df_joined$n_pos / df_joined$n_edges
df_all <- df_joined %>%
  group_by(dataset, partition, selection, positive_lim) %>%
  summarise(
    n_pos = sum(n_pos),
    .groups = "drop"
  ) %>%
  mutate(
    set_type = "total",
  ) %>%
  left_join(df_edge_count_dir, by = c("dataset", "positive_lim"))

df_all[, "percent_retained"] <- 100 * df_all$n_pos / df_all$n_edges

df_plot <- bind_rows(df_joined, df_all)
df_plot$set_type <- factor(df_plot$set_type, levels = c("total", "train", "val", "test"))

g <- ggplot(df_plot, aes(x = partition, y = selection)) +
  geom_tile(
    aes(
      fill = percent_retained
    )
  ) +
  geom_text(
    aes(
      label = sprintf("%.1f", percent_retained)
    ),
    color = "white", size = 3
  ) +
  scale_fill_gradient(high = "darkorange", low = "blue") +
  labs(
    title = "Positive edges retained in after partitioning",
    x = "Gene partition", y = "Selection", fill = " % edges retained"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -45), legend.position = "right") +
  facet_wrap(~ dataset + set_type, labeller = labeller(
    dataset = c(
      flat = "Combined",
      ms = "MS",
      y2h = "Y2H"
    ),
    set_type = c(
      total = "Total",
      train = "Train",
      val = "Validation",
      test = "Test"
    )
  )) +
  scale_x_discrete(labels = function(x) {
    case_when(
      str_remove(x, "\\..*") == "maxpos" ~ "MaxPos",
      str_remove(x, "\\..*") == "sequencesimilarity" ~ "SeqSim",
      TRUE ~ x
    )
  })


ggsave("percent_retained.png", g, height = 6, width = 6)



df_prob_dist <- read.table(
  "prediction_dist_ms_loose_directional_ss.csv",
  sep = "\t", header = TRUE
)

g <- ggplot(
  df_prob_dist,
  aes(x = value, fill = negative_data)
) +
  geom_density(alpha = .8) +
  theme_bw() +
  scale_fill_manual(
    values = c("darkorange", "blue"),
    labels = c(
      "HCNI", "Non-observed"
    )
  ) +
  labs(
    title = "Predicted probabilities HCNI vs non-observed\nMS dataset, loose selection, bait-prey network",
    x = "Prediction probability",
    fill = "Negative data"
  ) + theme(
    legend.position = "bottom"
  )

ggsave("prediction_probabilities.png", g, height = 4, width = 5)
