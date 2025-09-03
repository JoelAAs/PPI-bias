library(ggplot2)
library(tidyverse)
library(magrittr)

## Input
df_a_lesser = read.table(
  "work_folder/analysis/GO/abundance_jaccard_lesser.csv",
  sep="\t",
  header=T
)
df_a_lesser$model <- "abundance_aware"
df_a_greater = read.table(
  "work_folder/analysis/GO/abundance_jaccard_greater.csv",
  sep="\t",
  header=T
)
df_a_greater$model <- "abundance_aware"

df_flat_less= "work_folder/analysis/GO/flat_jaccard_lesser.csv"
df_flat_less %>%  read.table(sep="\t", header=T) -> less_flat
less_flat$limit <- "upper_bound_pod"
less_flat$model <- "flat"


df_flat_greater = "work_folder/analysis/GO/flat_jaccard_greater.csv"
df_flat_greater %>%  read.table(sep="\t", header=T) -> greater_flat
greater_flat$limit <- "lower_bound_pod"
greater_flat$model <- "flat"



df_expected_observed = bind_rows(
  df_a_lesser,
  df_a_greater,
  df_flat_less,
  df_flat_greater
)

go_jaccard <- ggplot(
  df_expected_observed,
  aes(
    x=value
  )
) +
  geom_point(aes(y = cc_ji_avg, color = "CC Jaccard")) +
  geom_point(aes(y = bp_ji_avg, color = "BP Jaccard")) +
  geom_point(aes(y = mf_ji_avg, color = "MF Jaccard")) +
  facet_wrap(model ~ limit,
             labeller = labeller(
               limit =
                 c("lower_bound_pod" = "Pmin > POD",
                   "upper_bound_pod" = "Pmax < POD")
             )) 
  labs(
    x = "Probability of detection",
    y = "Average Jaccard index",
    title = "GO jaccard index vs POD",
    color="GO category"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"))

ggsave("work_folder/plots/GO/jaccard_GO_vs_POD.png",
       go_jaccard,
       dpi=300,
       height=7,
       width=7
)



prob_intersect <- ggplot(
  df_expected_observed,
  aes(
    x=value
  )
) +
  geom_point(aes(y = cc_intersect_avg, color = "CC intersect")) +
  geom_point(aes(y = bf_intersect_avg, color = "BP intersect")) +
  geom_point(aes(y = mf_intersect_avg, color = "MF intersect")) +
  facet_wrap(model ~ limit,
             labeller = labeller(
               limit =
                 c("lower_bound_pod" = "Pmin > POD",
                   "upper_bound_pod" = "Pmax < POD")
             )) 
labs(
  x = "Probability of detection",
  y = "Average Jaccard index",
  title = "GO jaccard index vs POD",
  color="GO category"
) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"))

ggsave("work_folder/plots/GO/intersect_GO_vs_POD.png",
       prob_intersect,
       dpi=300,
       height=7,
       width=7
)

