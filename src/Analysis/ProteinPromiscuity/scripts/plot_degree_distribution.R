library(ggplot2)
library(dplyr)
library(tidyr)

read_detectability <- function(path) {
  df <- read.table(path, sep = "\t", header = TRUE)
  names(df)[1] <- "id"
  df
}

df <- inner_join(
  read_detectability(snakemake@input$ms_bait_detectability),
  read_detectability(snakemake@input$ms_prey_detectability),
  by = "id",
  suffix = c("_bait", "_prey")
) %>%
  mutate(detectability_mean = (detectability_bait + detectability_prey) / 2)

df_long <- df %>%
  select(id, bait = detectability_bait, prey = detectability_prey,
         mean = detectability_mean) %>%
  pivot_longer(-id, names_to = "role", values_to = "detectability") %>%
  mutate(role = factor(role, levels = c("bait", "prey", "mean")))

g <- ggplot(df_long, aes(x = detectability, fill = role)) +
  geom_histogram(
    bins = 30, position = position_dodge(preserve = "single"),
    colour = "black", linewidth = 0.2
  ) +
  scale_x_log10() +
  scale_fill_manual(
    values = c(bait = "forestgreen", prey = "darkorange", mean = "grey50")
  ) +
  labs(
    x = "log(MS detectability)",
    y = "Number of proteins",
    fill = NULL
  ) +
  theme_bw()

ggsave(snakemake@output$degree_distribution_plot, g,
       dpi = 300, height = 4, width = 6)
