library(ggplot2)
library(dplyr)

read_partner_effect <- function(path) {
  df <- read.table(path, sep = "\t", header = TRUE)
  names(df)[1] <- "id"
  # p = detected fraction of the tests the protein takes part in,
  # n_pairs = distinct partners it was tested against
  df[, c("id", "p", "n_pairs")]
}

read_dataset <- function(bait_path, prey_path, dataset) {
  bait <- read_partner_effect(bait_path)
  prey <- read_partner_effect(prey_path)

  cat(sprintf(
    "%s: keeping %d protein(s) tested in both roles, dropping %d bait-only and %d prey-only\n",
    dataset, length(intersect(bait$id, prey$id)),
    length(setdiff(bait$id, prey$id)), length(setdiff(prey$id, bait$id))
  ))

  inner_join(bait, prey, by = "id", suffix = c("_bait", "_prey")) %>%
    mutate(dataset = dataset)
}

df <- bind_rows(
  read_dataset(snakemake@input$ms_bait, snakemake@input$ms_prey, "ms"),
  read_dataset(snakemake@input$y2h_bait, snakemake@input$y2h_prey, "y2h")
)

stats <- df %>%
  group_by(dataset) %>%
  summarise(
    rho = cor(p_bait, p_prey, method = "spearman"),
    # ties are expected, exact p is not attempted
    rho_p = suppressWarnings(
      cor.test(p_bait, p_prey, method = "spearman", exact = FALSE)$p.value),
    # paired on the log2 scale, so the test matches the reported fold change
    log2_fc = median(log2(p_bait / p_prey)),
    wilcox_p = suppressWarnings(
      wilcox.test(log2(p_bait), log2(p_prey), paired = TRUE)$p.value),
    n = n(),
    .groups = "drop"
  )

stat_labels <- stats %>%
  mutate(label = sprintf(
    "Spearman rho = %.3f (p = %.2e)\nlog2 FC (bait/prey) = %.3f\nWilcoxon p = %.2e (n = %d)",
    rho, rho_p, log2_fc, wilcox_p, n
  ))

lims <- range(c(df$p_bait, df$p_prey), na.rm = TRUE)

g <- ggplot(df, aes(x = p_bait, y = p_prey)) +
  geom_point(aes(color = dataset), alpha = 0.4, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_label(
    data = stat_labels,
    aes(x = 0, y = Inf, label = label),
    hjust = -0.05, vjust = 1.2, inherit.aes = FALSE, size = 3,
    label.size = 0.3, fill = "white"
  ) +
  facet_wrap(~dataset, ncol = 2, labeller = as_labeller(toupper)) +
  scale_x_log10(limits = lims) +
  scale_y_log10(limits = lims) +
  scale_color_manual(values = c(ms = "forestgreen", y2h = "darkorange"), guide = "none") +
  labs(x = "log(Detectability as bait)", y = "log(Detectability as prey)") +
  theme_bw()

ggsave(snakemake@output$detectability, g, dpi = 300, height = 4, width = 8)

# Detectability on a common scale: the role bias (median log2 FC, bait over prey) is
# divided out of the bait side, putting it on the prey scale, and the two are averaged.
adjusted <- df %>%
  left_join(stats[, c("dataset", "log2_fc")], by = "dataset") %>%
  mutate(
    p_bait_adjusted = p_bait * 2 ^ -log2_fc,
    adjusted_detectability = (p_bait_adjusted + p_prey) / 2,
    min_n_pairs = pmin(n_pairs_bait, n_pairs_prey)  # partners backing the weaker role
  )

write_adjusted <- function(ds, path) {
  adjusted %>%
    filter(dataset == ds) %>%
    select(id, p_bait, p_prey, p_bait_adjusted, adjusted_detectability, min_n_pairs) %>%
    arrange(desc(adjusted_detectability)) %>%
    write.table(path, sep = "\t", row.names = FALSE, quote = FALSE)
}

write_adjusted("ms", snakemake@output$ms_adjusted_detectability)
write_adjusted("y2h", snakemake@output$y2h_adjusted_detectability)
