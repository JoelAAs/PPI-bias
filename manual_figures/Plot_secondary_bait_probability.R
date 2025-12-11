library(ggplot2)
library(tidyverse)
library(reshape2)
library(stringr)


df_entropy_ms = read.table("work_folder/per_gene/analysis/negatome/test_entropy_ms_limit_3.csv", sep="\t", header=T)
df_entropy_y2h = read.table("work_folder/per_gene/analysis/negatome/test_entropy_y2h_limit_3.csv", sep="\t", header=T)
df_entropy_ms$method <- "AP/IP - MS"
df_entropy_y2h$method <- "Y2H"

df_entropy <- bind_rows(df_entropy_ms, df_entropy_y2h)

entropy_plot <- ggplot(df_entropy,
  aes(
    x=as.factor(n_tested),
    y=2^pair_entropy- n_tested,
    fill=method
  )) +
  geom_boxplot() +
  theme_bw() +
  labs(
    x="Times protein pair tested",
    y=expression(paste(2^H, "- total tests")),
    fill="Detection method"
  ) +
  scale_fill_manual(
    values = c("darkorange", "blue")
  ) +
  theme(
    legend.position = "bottom"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed")

max_prob_plot <- ggplot(df_entropy,
                       aes(
                         x=as.factor(n_tested),
                         y=max_bait_probability - 1/n_tested,
                         fill=method
                       )) +
  geom_boxplot() +
  theme_bw() +
  labs(
    x="Times protein pair tested",
    y=expression(paste(P(b[i])[max], "-",  1/N[Tests])),
    fill="Detection method"
  ) +
  scale_fill_manual(
    values = c("darkorange", "blue")
  ) +
  theme(
    legend.position = "bottom"
  ) +
  geom_hline(yintercept = 0, linetype="dashed")


ggsave("work_folder/per_gene/plotting/negatome/pair_testing_entropy.png",
       entropy_plot,
       dpi=300,
       height=3,
       width=5)

ggsave("work_folder/per_gene/plotting/negatome/pair_testing_probability.png",
       max_prob_plot,
       dpi=300,
       height=3,
       width=5)

max_prob_plot <- df_entropy %>%  group_by(
  method, n_tested
) %>% 
  summarise(
    mean_entropy = mean(pair_entropy),
    mean_max_p = mean(max_bait_probability),
    number_of_rows = n()
  ) %>% ungroup()
