library(ggplot2)
library(dplyr)

statistic_label <- function(x) {
  case_when(
    x == "oe_negative" ~ "Shared negative",
    x == "oe_positive" ~ "Shared positive",
    TRUE ~ x
  )
}

read_one <- function(path) {
  read.table(path, sep = "\t", header = TRUE)
}

df <- bind_rows(lapply(snakemake@input, read_one)) %>%
  filter(statistic %in% c("oe_negative", "oe_positive"), estimate_type == "raw") %>%
  mutate(
    statistic_label = factor(
      statistic_label(statistic),
      levels = c("Shared positive", "Shared negative")
    )
  )

g <- ggplot(df, aes(x = value, y = statistic_label)) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(size = 2) +
  labs(
    title = "MS/Y2H concordance (any interaction observed)",
    x = "O/E ratio",
    y = ""
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 9)
  )

ggsave(snakemake@output[[1]], g, dpi = 300, height = 3, width = 6)
