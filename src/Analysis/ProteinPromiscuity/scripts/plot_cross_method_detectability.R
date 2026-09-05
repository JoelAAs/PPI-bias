library(ggplot2)
library(dplyr)

capitalise <- function(x) paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))

read_detectability <- function(path) {
  df <- read.table(path, sep = "\t", header = TRUE)
  names(df)[1] <- "id"
  df
}

join_role <- function(ms_path, y2h_path, role) {
  inner_join(
    read_detectability(ms_path),
    read_detectability(y2h_path),
    by = "id",
    suffix = c("_ms", "_y2h")
  ) %>%
    mutate(role = role)
}

df <- bind_rows(
  join_role(snakemake@input$ms_bait_detectability, snakemake@input$y2h_bait_detectability, "bait"),
  join_role(snakemake@input$ms_prey_detectability, snakemake@input$y2h_prey_detectability, "prey")
)

stat_labels <- df %>%
  group_by(role) %>%
  summarise(
    rho = cor(detectability_ms, detectability_y2h, method = "spearman"),
    log2_enrichment = median(log2(detectability_ms / detectability_y2h)),
    wilcox_p = wilcox.test(detectability_ms, detectability_y2h, paired = TRUE)$p.value,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf(
    "Spearman rho = %.3f\nlog2 FC (ms/y2h) = %.3f\nWilcoxon p = %.2e (n = %d)",
    rho, log2_enrichment, wilcox_p, n
  ))

lims_ms_y2h <- range(c(df$detectability_ms, df$detectability_y2h), na.rm = TRUE)

g <- ggplot(df, aes(x = detectability_ms, y = detectability_y2h)) +
  geom_point(alpha = 0.4, size = 1, color="blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_label(
    data = stat_labels,
    aes(x = 0, y = Inf, label = label),
    hjust = -0.05, vjust = 1.2, inherit.aes = FALSE, size = 3,
    linewidth = 0.3, fill = "white"
  ) +
  facet_wrap(~role, ncol = 2, labeller = as_labeller(capitalise)) +
  scale_x_log10(limits = lims_ms_y2h) +
  scale_y_log10(limits = lims_ms_y2h) +
  labs(x = "log(MS detectability)", y = "log(Y2H detectability)") +
  theme_bw() +
  theme(
  )

ggsave(snakemake@output$ms_y2h_plot, g, dpi = 300, height = 4, width = 8)

join_dataset <- function(bait_path, prey_path, dataset) {
  inner_join(
    read_detectability(bait_path),
    read_detectability(prey_path),
    by = "id",
    suffix = c("_bait", "_prey")
  ) %>%
    mutate(dataset = dataset)
}

df_bp <- bind_rows(
  join_dataset(snakemake@input$ms_bait_detectability, snakemake@input$ms_prey_detectability, "ms"),
  join_dataset(snakemake@input$y2h_bait_detectability, snakemake@input$y2h_prey_detectability, "y2h")
)

stat_labels_bp <- df_bp %>%
  group_by(dataset) %>%
  summarise(
    rho = cor(detectability_bait, detectability_prey, method = "spearman"),
    # axes are log10, so Pearson is taken on the logged values
    pearson_r = cor(log10(detectability_bait), log10(detectability_prey)),
    log2_enrichment = median(log2(detectability_bait / detectability_prey)),
    wilcox_p = wilcox.test(detectability_bait, detectability_prey, paired = TRUE)$p.value,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf(
    paste0(
      "Spearman rho = %.3f\nPearson r (log10) = %.3f\n",
      "log2 FC (bait/prey) = %.3f\nWilcoxon p = %.2e (n = %d)"
    ),
    rho, pearson_r, log2_enrichment, wilcox_p, n
  ))

lims_bait_prey <- range(c(df_bp$detectability_bait, df_bp$detectability_prey), na.rm = TRUE)

g_bp <- ggplot(df_bp, aes(x = detectability_bait, y = detectability_prey)) +
  geom_point(aes(color = dataset), alpha = 0.4, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_label(
    data = stat_labels_bp,
    aes(x = 0, y = 0, label = label),
    hjust = -0.05, vjust = -0.2, inherit.aes = FALSE, size = 3,
    linewidth = 0.3, fill = "white"
  ) +
  facet_wrap(~dataset, ncol = 2, labeller = as_labeller(toupper)) +
  scale_x_log10(limits = lims_bait_prey) +
  scale_y_log10(limits = lims_bait_prey) +
  scale_color_manual(values = c(ms = "forestgreen", y2h = "darkorange"), guide = "none") +
  labs(x = "log(Bait detectability)", y = "log(Prey detectability)") +
  theme_bw()

ggsave(snakemake@output$bait_prey_plot, g_bp, dpi = 300, height = 4, width = 8)
