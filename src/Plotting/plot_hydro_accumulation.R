library(ggplot2)
library(tidyverse)
library(magrittr)


min_pairs = 200

## Input
df_a_lesser = read.table(
  "work_folder/analysis/hydrophobicity/abundance_netsurfp_lesser.csv",
  sep="\t",
  header=T
)
df_a_lesser$model <- "abundance_aware"


df_a_greater = read.table(
  "work_folder/analysis/hydrophobicity/abundance_netsurfp_greater.csv",
  sep="\t",
  header=T
)
df_a_greater$model <- "abundance_aware"

df_flat_less= "work_folder/analysis/hydrophobicity/flat_netsurfp_lesser.csv"
df_flat_less %>%  read.table(sep="\t", header=T) -> less_flat
less_flat$limit <- "upper_bound_pod"
less_flat$model <- "flat"


df_flat_greater = "work_folder/analysis/hydrophobicity/flat_netsurfp_greater.csv"
df_flat_greater %>%  read.table(sep="\t", header=T) -> greater_flat
greater_flat$limit <- "lower_bound_pod"
greater_flat$model <- "flat"
prob  <- function(x) 1/(1+exp(-x))
logit <- function(x) -log(1/x-1)


abundance_df = bind_rows(
  df_a_lesser,
  df_a_greater
) %>% filter(number_of_pairs > min_pairs)

abundance_df$thsa_delta_avg <- abundance_df$thsa_delta_avg/max(abundance_df$thsa_delta_avg)
abundance_df$tasa_delta_avg <- abundance_df$tasa_delta_avg/max(abundance_df$tasa_delta_avg)
abundance_df$rhsa_delta_avg <- abundance_df$rhsa_delta_avg/max(abundance_df$rhsa_delta_avg)

flat_df = bind_rows(
  less_flat,
  greater_flat
) %>% filter(number_of_pairs > min_pairs)

flat_df$thsa_delta_avg <- flat_df$thsa_delta_avg/max(flat_df$thsa_delta_avg)
flat_df$tasa_delta_avg <- flat_df$tasa_delta_avg/max(flat_df$tasa_delta_avg)
flat_df$rhsa_delta_avg <- flat_df$rhsa_delta_avg/max(flat_df$rhsa_delta_avg)

df_expected_observed = bind_rows(
  df_a_lesser,
  df_a_greater,
  less_flat,
  greater_flat
) %>% filter(number_of_pairs > 200)


### Abundance
go_rsha_a <- ggplot(
  abundance_df,
  aes(
    x=prob(value)
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
    title = "Hydrophobicity vs abundance POD",
    color="Hydrophobicity"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"))


### Flat
go_hydro_flat <- ggplot(
  flat_df,
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
    title = "Hydrophobicity vs flat POD",
    color="Hydrophobicity"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"))






add_label <- function(plot, label) {
  arrangeGrob(plot, top = textGrob(label, x = unit(0.05, "npc"), y = unit(0.9, "npc"),
                                   just = c("left", "top"),
                                   gp = gpar(fontsize = 16, fontface = "bold")))
}



flat_GO <- grid.arrange(
  add_label(go_rsha_a, "A"),
  add_label(go_hydro_flat, "B"),
  nrow=2)



ggsave("work_folder/plots/Hydro/hydro_accum.png",
       flat_GO,
       dpi=300,
       height=6,
       width=4
)
