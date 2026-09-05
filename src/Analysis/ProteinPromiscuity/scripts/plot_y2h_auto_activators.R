library(ggplot2)
library(dplyr)

capitalise <- function(x) paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))

read_keyed <- function(path) {
  df <- read.table(path, sep = "\t", header = TRUE)
  names(df)[1] <- "id"
  df
}

join_role <- function(dispersion_path, detectability_path) {
  inner_join(
    read_keyed(dispersion_path),
    read_keyed(detectability_path),
    by = "id"
  )
}

df <- bind_rows(
  join_role(snakemake@input$prey_auto_activators, snakemake@input$y2h_prey_detectability),
  join_role(snakemake@input$bait_auto_activators, snakemake@input$y2h_bait_detectability)
)

stat_labels <- df %>%
  group_by(role) %>%
  summarise(
    rho = cor(detectability, dispersion_ratio, method = "spearman"),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("Spearman rho = %.3f (n = %d)", rho, n))

# a protein detected at the same rate by every partner has obs_var == 0; floor
# those onto the log axis rather than dropping them, they are the sticky ones
floor_ratio <- min(df$dispersion_ratio[df$dispersion_ratio > 0]) / 2
df$dispersion_ratio_floored <- pmax(df$dispersion_ratio, floor_ratio)

g <- ggplot(df, aes(y = detectability, x = dispersion_ratio_floored)) +
  geom_point(alpha = 0.4, size = 1, color = "darkorange") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_label(
    data = stat_labels,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.05, vjust = 1.2, inherit.aes = FALSE, size = 3,
    linewidth = 0.3, fill = "white"
  ) +
  facet_wrap(~role, ncol = 2, labeller = as_labeller(capitalise)) +
  scale_y_log10() +
  labs(
    y = "log(Y2H detectability)",
    x = "Dispersion ratio: Var(observed) / Var(expected under uniform rate)"
  ) +
  theme_bw()

ggsave(snakemake@output$auto_activator_plot, g, dpi = 300, height = 4, width = 8)
