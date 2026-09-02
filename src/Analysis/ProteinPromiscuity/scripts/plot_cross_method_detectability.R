library(ggplot2)
library(dplyr)

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

g <- ggplot(df, aes(x = detectability_ms, y = detectability_y2h)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~role, ncol = 2) +
  labs(x = "MS detectability", y = "Y2H detectability") +
  theme_bw()

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

g_bp <- ggplot(df_bp, aes(x = detectability_bait, y = detectability_prey)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~dataset, ncol = 2) +
  labs(x = "Bait detectability", y = "Prey detectability") +
  theme_bw()

ggsave(snakemake@output$bait_prey_plot, g_bp, dpi = 300, height = 4, width = 8)
