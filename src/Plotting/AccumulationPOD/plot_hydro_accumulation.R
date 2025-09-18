library(ggplot2)
library(tidyverse)
library(magrittr)

prob  <- function(x) 1/(1+exp(-x))
logit <- function(x) -log(1/x-1)

### args
args <- commandArgs(trailingOnly = TRUE)

greater_hydro_accumulation <- args[1]
lesser_hydro_accumulation  <- args[2]
name                       <- args[3]
output                     <- args[4]

### Input
df_greater = read.table(
  greater_hydro_accumulation,
  sep="\t",
  header=T
)

df_lesser = read.table(
  lesser_hydro_accumulation,
  sep="\t",
  header=T
)

hydro_df = bind_rows(
  df_greater,
  df_lesser
)

hydro_df$thsa_netsurfp2_delta <- hydro_df$sum_thsa_netsurfp2_delta/hydro_df$non_na_pairs_thsa_netsurfp2_delta
hydro_df$tasa_netsurfp2_delta <- hydro_df$sum_tasa_netsurfp2_delta/hydro_df$non_na_pairs_tasa_netsurfp2_delta
hydro_df$rhsa_netsurfp2_delta <- hydro_df$sum_rhsa_netsurfp2_delta/hydro_df$non_na_pairs_rhsa_netsurfp2_delta

### plot
plot_rhsa <- ggplot(
  hydro_df,
  aes(
    x=value
  )
) +
  geom_point(aes(y = thsa_delta_avg, color = "THSA")) +
  geom_point(aes(y = tasa_delta_avg, color = "TASA")) +
  geom_point(aes(y = rhsa_delta_avg, color = "RHSA")) +
  geom_line(aes(y = thsa_delta_avg, color = "THSA")) +
  geom_line(aes(y = tasa_delta_avg, color = "TASA")) +
  geom_line(aes(y = rhsa_delta_avg, color = "RHSA")) +
  scale_color_discrete(labels = c(
    expression(Delta ~ THSA),
    expression(Delta ~ TASA),
    expression(Delta ~ RHSA)
    )) +
  facet_wrap(. ~ limit,
             labeller = labeller(
               limit =
                 c("lower_bound_pod" = "mean(Pmin > POD)",
                   "upper_bound_pod" = "mean(Pmax < POD)")
             ), scales = "free_x") +
  labs(
    x = "Probability of detection",
    y = "Normalised difference",
    title = paste("Hydrophobicity vs POD:", name),
    color="Hydrophobicity"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"))



ggsave(output,
       plot_rhsa,
       dpi=300,
       height=6,
       width=4
)
