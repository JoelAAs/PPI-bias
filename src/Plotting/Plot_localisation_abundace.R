library(ggplot2)
library(tidyverse)
library(magrittr)

## Input
df_neg = read.table(
  "work_folder/analysis/abundance_aware/localisation/probability_match_abundance_less.csv",
  sep="\t",
  header=T
  )
df_neg$model <- "abundance_aware"
df_hci = read.table(
  "work_folder/analysis/abundance_aware/localisation/probability_match_abundance_greater.csv",
  sep="\t",
  header=T
)
df_hci$model <- "abundance_aware"

localisation_match_greater = "work_folder/inferred_search_space/analysis/bias_reduced_ppis/localisation_p_estimated_protein_pairs_greater.csv"
localisation_match_greater %>%  read.table(sep="\t", header=T) -> greater_flat
greater_flat$limit <- "lower_bound_pod"
greater_flat$model <- "flat"

localisation_match_less= "work_folder/inferred_search_space/analysis/bias_reduced_ppis/localisation_p_estimated_protein_pairs_less.csv"
localisation_match_less %>%  read.table(sep="\t", header=T) -> less_flat
less_flat$limit <- "upper_bound_pod"
less_flat$model <- "flat"


df_expected_observed = bind_rows(
  df_neg,
  df_hci,
  greater_flat,
  less_flat
)


df_expected_observed$percent_change <- (df_expected_observed$observed - df_expected_observed$expected)/df_expected_observed$expected
prob <- function(x) 1/(1 + exp(-x))
df_expected_observed$prob <- sapply(df_expected_observed$value, prob)



prob_localisation <- ggplot(
  df_expected_observed,
  aes(
    x=prob,
    y=percent_change,
    color=model
  )
) +
  geom_point() +
  geom_line() + 
  facet_wrap(. ~ limit,
             labeller = labeller(
               limit =
                 c("lower_bound_pod" = "P(m| Pmin > POD)",
                   "upper_bound_pod" = "P(m| Pmax < POD)")
  )) +
  labs(
    x = "Probability of detection",
    y = "% change in expected localisation match",
    title = "Expected localisation matches divegrence"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"))

ggsave("work_folder/plots/localisation/POD_vs_localisation_match.png",
       prob_localisation,
       dpi=300,
       height=4,
       width=7
       )

ggplot(
  df_expected_observed,
  aes(
    x=value,
    y=number_of_pairs,
    color=model
  )
) +
  geom_point() +
  geom_line() + 
  facet_wrap(. ~ limit,
             labeller = labeller(
               limit =
                 c("lower_bound_pod" = "P(m| Pmin > POD)",
                   "upper_bound_pod" = "P(m| Pmax < POD)")
  )) +
  labs(
    x = "Probability of detection",
    y = "% change in expected localisation match",
    title = "Expected localisation matches divegrence"
  ) +
  theme_bw() +
  theme(legend.position = "right",
        strip.text = element_text(size = 10, face = "bold"))
