library(ggplot2)
library(tidyverse)
library(magrittr)

### args
args <- commandArgs(trailingOnly = TRUE)

greater_colocalisation_accumulation <- args[1]
lesser_colocalisation_accumulation  <- args[2]
name                       <- args[3]
output                     <- args[4]

### Input
df_greater = read.table(
  greater_colocalisation_accumulation,
  sep="\t",
  header=T
)
df_lesser = read.table(
  lesser_colocalisation_accumulation,
  sep="\t",
  header=T
)

colocalisation_df = bind_rows(
  df_lesser,
  df_greater
)
prob  <- function(x) 1/(1+exp(-x))
logit <- function(x) -log(1/x-1)
if (name=="abundance_mcmc") {
    colocalisation_df$value <- prob(colocalisation_df$value)
}



colocalisation_df$percent_change <- (colocalisation_df$sum_localisation_match - colocalisation_df$sum_match_probability)/colocalisation_df$sum_match_probability

prob_localisation <- ggplot(
  colocalisation_df,
  aes(
    x=value,
    y=percent_change,
  )
) +
  geom_point(color="red") +
  geom_line() + 
  facet_wrap(. ~ limit,
             labeller = labeller(
               limit =
                 c("lower_bound_pod" = "P(m| Pmin > POD)",
                   "upper_bound_pod" = "P(m| Pmax < POD)")
  ), scales="free_x") +
  labs(
    x = "Probability of detection",
    y = "% different in expected localisation match",
    title = paste("Observed vs Expected colocalisation:", name)
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"))

ggsave(output,
       prob_localisation,
       dpi=300,
       height=4,
       width=7
       )
