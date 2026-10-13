library(ggplot2)
library(dplyr)

bp_df <- read.table(snakemake@input$bp_df, sep = "\t", header = TRUE)

bp_df %>%
  filter(detection_method %in% snakemake@params$ms_methods) %>%
  mutate(
    uniprot_id_bait = sapply(
      uniprot_id_bait, function(x) strsplit(x, "-")[[1]][1]),
    uniprot_id_prey = sapply(
      uniprot_id_prey, function(x) strsplit(x, "-")[[1]][1])
  ) %>%
  filter(uniprot_id_bait != uniprot_id_prey) -> bp_df

study_counts <- bp_df %>%
  mutate(study = paste(pubmed_id, detection_method, sep = "_")) %>%
  group_by(study) %>%
  summarise(
    n_bait = n_distinct(uniprot_id_bait),
    n_prey = n_distinct(uniprot_id_prey),
    .groups = "drop"
  )

# Michaelis-Menten: y = Vm * x / (K + x); Vm is the saturation (asymptote) of prey
mm_fit <- nls(n_prey ~ SSmicmen(n_bait, Vm, K), data = study_counts)
coefs <- summary(mm_fit)$coefficients
rss <- sum(residuals(mm_fit)^2)
tss <- sum((study_counts$n_prey - mean(study_counts$n_prey))^2)
vm <- coefs["Vm", "Estimate"]
vm_se <- coefs["Vm", "Std. Error"]
fit_stats <- data.frame(
  n_studies = nrow(study_counts),
  Vm_saturation = vm,
  Vm_ci_low = vm - qt(0.975, df.residual(mm_fit)) * vm_se,
  Vm_ci_high = vm + qt(0.975, df.residual(mm_fit)) * vm_se,
  K_half_saturation_baits = coefs["K", "Estimate"],
  K_se = coefs["K", "Std. Error"],
  pseudo_R2 = 1 - rss / tss,
  RMSE = sqrt(rss / nrow(study_counts)),
  AIC = AIC(mm_fit),
  max_observed_prey = max(study_counts$n_prey),
  frac_of_saturation_max_observed = max(study_counts$n_prey) / vm
)
write.table(fit_stats, snakemake@output$fit_stats, sep = "\t",
            quote = FALSE, row.names = FALSE)

g <- ggplot(study_counts, aes(x = n_bait, y = n_prey)) +
  geom_point(alpha = 0.4) +
  geom_smooth(
    method = "nls",
    formula = y ~ SSmicmen(x, Vm, K),  # Michaelis-Menten: y = Vm * x / (K + x)
    se = FALSE,
    color = "steelblue"
  ) +
  geom_hline(yintercept = vm, linetype = "dashed", color = "firebrick") +
  annotate(
    "text", x = 0, y = vm, hjust = 0, vjust = -0.5, color = "firebrick",
    label = sprintf("Saturation Vm = %.0f prey (pseudo-R2 = %.2f)",
                    vm, fit_stats$pseudo_R2)
  ) +
  labs(x = "Unique baits", y = "Unique prey") +
  theme_bw()

ggsave(snakemake@output$plot, g, dpi = 300, height = 5, width = 6)
