library(ggplot2)
library(tidyverse)
library(magrittr)
library(ggExtra)

df = read.table(
  "work_folder/inferred_search_space/analysis/cell_line/bait_prior_long.csv",
  sep="\t",
  header=T
)

df_or <- read.table(
  "work_folder/inferred_search_space/analysis/cell_line/bait_wise_prey_plotting.csv",
  sep="\t",
  header=T
)

df <- df %>% 
  filter(n_tested > 4) %>% 
  filter(n_tested != total_tested)

df_or_p <- merge(
  df,
  df_or,
  by=c("gene_name_prey", "cl_id")
  )

df$delta_p = df$prey_p - df$total_observed/df$total_tested

ggplot(
  df,
  aes(
    x=n_tested,
    y=n_observed,
    colour = cl_id
  )
) + geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  facet_wrap(. ~cl_id)

ggplot(
  df,
  aes(
    x=delta_p,
    y=n_tested,
    color=cl_id
  )
) + geom_point() +
facet_wrap(. ~cl_id)


ggplot(
  df,
    aes(
    x = total_tested,
    y = total_observed
  )
) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_log10() +
  scale_y_log10() 
  

p = ggplot(
  df_or_p,
  aes(
    x = delta_p,
    y = log(odds_ratio),
    color=cl_id
  )
) + geom_point() +
  facet_wrap(cl_id ~.) +
  theme_bw() +
  theme(
    legend.position="bottom"
    )

ggsave("p_vs_or.png", p, dpi=300, height=4, width=4*1.618)


ggplot(
  df,
  aes(x=n_tested,
      fill=cl_id)
) + geom_histogram(position="dodge", bins=40) +
  theme_bw()
