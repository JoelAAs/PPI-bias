library(ggplot2)
library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)

min_pairs = 200

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
prob  <- function(x) 1/(1+exp(-x))
logit <- function(x) -log(1/x-1)


abundance_df = bind_rows(
  df_a_lesser,
  df_a_greater
) %>% filter(number_of_pairs > min_pairs)

flat_df = bind_rows(
  less_flat,
  greater_flat
) %>% filter(number_of_pairs > min_pairs)

df_expected_observed = bind_rows(
  df_a_lesser,
  df_a_greater,
  less_flat,
  greater_flat
) %>% filter(number_of_pairs > 200)


### Abundance
go_jaccard_a <- ggplot(
  abundance_df,
  aes(
    x=prob(value)
  )
) +
  geom_point(aes(y = cc_ji_avg, color = "CC Jaccard")) +
  geom_point(aes(y = bp_ji_avg, color = "BP Jaccard")) +
  geom_point(aes(y = mf_ji_avg, color = "MF Jaccard")) +
  facet_wrap(. ~ limit,
             labeller = labeller(
               limit =
                 c("lower_bound_pod" = "mean(Pmin > POD)",
                   "upper_bound_pod" = "mean(Pmax < POD)")
             ), scales = "free_x") +
  labs(
    x = "Probability of detection",
    y = "Average Jaccard index",
    title = "GO jaccard index vs abundance POD",
    color="GO category"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"))


prob_intersect_a <- ggplot(
  abundance_df,
  aes(
    x=prob(value)
  )
) +
  geom_point(aes(y = cc_intersect_avg, color = "CC intersect")) +
  geom_point(aes(y = bp_intersect_avg, color = "BP intersect")) +
  geom_point(aes(y = mf_intersect_avg, color = "MF intersect")) +
  facet_wrap(. ~ limit,
             labeller = labeller(
               limit =
                 c("lower_bound_pod" = "Pmin > POD",
                   "upper_bound_pod" = "Pmax < POD")
             ), scales = "free_x") +
labs(
  x = "Probability of detection",
  y = "Avergage length intersection",
  title = "GO Intersection index vs abundance POD",
  color="GO category"
) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"))


add_label <- function(plot, label) {
  arrangeGrob(plot, top = textGrob(label, x = unit(0.05, "npc"), y = unit(0.9, "npc"),
                                   just = c("left", "top"),
                                   gp = gpar(fontsize = 16, fontface = "bold")))
}

abundance_GO <- grid.arrange(
  add_label(go_jaccard_a, "A"),
  add_label(prob_intersect_a, "B"),
  nrow=2)

ggsave("work_folder/plots/GO/Abundance_GO.png",
       abundance_GO,
       dpi=300,
       height=6,
       width=6
)

### Flat
go_jaccard_flat <- ggplot(
  flat_df,
  aes(
    x=value
  )
) +
  geom_point(aes(y = cc_ji_avg, color = "CC Jaccard")) +
  geom_point(aes(y = bp_ji_avg, color = "BP Jaccard")) +
  geom_point(aes(y = mf_ji_avg, color = "MF Jaccard")) +
  facet_wrap(. ~ limit,
             labeller = labeller(
               limit =
                 c("lower_bound_pod" = "mean(Pmin > POD)",
                   "upper_bound_pod" = "mean(Pmax < POD)")
             ), scales = "free_x") +
  labs(
    x = "Probability of detection",
    y = "Average Jaccard index",
    title = "GO jaccard index vs flat POD",
    color="GO category"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"))


prob_intersect_flat <- ggplot(
  flat_df,
  aes(
    x=value
  )
) +
  geom_point(aes(y = cc_intersect_avg, color = "CC intersect")) +
  geom_point(aes(y = bp_intersect_avg, color = "BP intersect")) +
  geom_point(aes(y = mf_intersect_avg, color = "MF intersect")) +
  facet_wrap(. ~ limit,
             labeller = labeller(
               limit =
                 c("lower_bound_pod" = "Pmin > POD",
                   "upper_bound_pod" = "Pmax < POD")
             ), scales = "free_x") +
  labs(
    x = "Probability of detection",
    y = "Avergage length intersection",
    title = "GO Intersection index vs flat POD",
    color="GO category"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"))


add_label <- function(plot, label) {
  arrangeGrob(plot, top = textGrob(label, x = unit(0.05, "npc"), y = unit(0.9, "npc"),
                                   just = c("left", "top"),
                                   gp = gpar(fontsize = 16, fontface = "bold")))
}



low_ji <- ggplot(
  flat_df %>% filter(limit =="upper_bound_pod") %>%  filter(value < 0.25),
  aes(
    x=logit(value)
  )
) +
  geom_point(aes(y = cc_ji_avg, color = "CC intersect")) +
  geom_point(aes(y = bp_ji_avg, color = "BP intersect")) +
  geom_point(aes(y = mf_ji_avg, color = "MF intersect")) +
  labs(
    x = "logit(Probability of detection)",
    y = "Average Jaccard index",
    color="GO category"
  ) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"))

low_inter <- ggplot(
  flat_df %>% filter(limit =="upper_bound_pod") %>%  filter(value < 0.25),
  aes(
    x=logit(value)
  )
) +
  geom_point(aes(y = cc_intersect_avg, color = "CC intersect")) +
  geom_point(aes(y = bp_intersect_avg, color = "BP intersect")) +
  geom_point(aes(y = mf_intersect_avg, color = "MF intersect")) +
  labs(
    x = "logit(Probability of detection)",
    y = "Avergage length intersection",
    color="GO category"
  ) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"))


add_label <- function(plot, label) {
  arrangeGrob(plot, top = textGrob(label, x = unit(0.05, "npc"), y = unit(0.9, "npc"),
                                   just = c("left", "top"),
                                   gp = gpar(fontsize = 16, fontface = "bold")))
}



flat_GO <- grid.arrange(
  add_label(go_jaccard_flat, "A"),
  add_label(prob_intersect_flat, "B"),
  add_label(grid.arrange(low_ji,low_inter, nrow=1), "C"),
  nrow=3)



ggsave("work_folder/plots/GO/flat_GO.png",
       flat_GO,
       dpi=300,
       height=10,
       width=6
)
