library(ggplot2)
library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)

min_pairs = 200

## Input
args <- commandArgs(trailingOnly = TRUE)

greater_do_accumulation <- args[1]
lesser_do_accumulation  <- args[2]
name                    <- args[3]
output_jacccard         <- args[4]
output_accumulation     <- args[5]

prob  <- function(x) 1/(1+exp(-x))
logit <- function(x) -log(1/x-1)

df_greater = read.table(
  greater_do_accumulation,
  sep="\t",
  header=T
)
df_lesser = read.table(
  lesser_do_accumulation,
  sep="\t",
  header=T
)

do_df = bind_rows(
  df_greater,
  df_lesser
)

if (name=="abundance_mcmc") {
    do_df$value <- prob(do_df$value)
    c_xlab = "Probability of detection"
} else {
    do_df[do_df$limit =="upper_bound_pod", "value"] = logit(do_df[do_df$limit =="upper_bound_pod", "value"])
    c_xlab = "Probability of detection -- logit(POD)"
}




do_df$ji <- do_df$sum_ji_do/do_df$non_na_pairs_ji_do
do_df$intersect_do <- do_df$sum_intersect_do/do_df$non_na_pairs_intersect_do

do_df$n_do_bait <- do_df$sum_n_do_bait/do_df$non_na_pairs_n_do_bait
do_df$n_do_prey <- do_df$sum_n_do_prey/do_df$non_na_pairs_n_do_prey

do_df$limit <- factor(do_df$limit, levels = c("upper_bound_pod", "lower_bound_pod"))

### do
plot_jaccard <- ggplot(
  do_df,
  aes(
    x=value,
    y=ji
  )
) +
  geom_point(aes(color = "DO intersect")) +
  geom_line(aes(color = "DO intersect")) +
  facet_wrap(. ~ limit,
             labeller = labeller(
               limit =
                 c("lower_bound_pod" = "mean(Pmin > POD)",
                   "upper_bound_pod" = "mean(Pmax < POD)")
             ), scales = "free_x") +
  labs(
    x = c_xlab,
    y = "Average Jaccard index",
    title = paste("do jaccard index vs POD:", name),
  ) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"))


prob_intersect <- ggplot(
  do_df,
  aes(
    x=value,
    y = intersect_do
  )
) +
  geom_point(aes(color = "DO intersect")) +
  geom_line(aes(color = "DO intersect")) +
  facet_wrap(. ~ limit,
             labeller = labeller(
               limit =
                 c("lower_bound_pod" = "Pmin > POD",
                   "upper_bound_pod" = "Pmax < POD")
             ), scales = "free_x") +
labs(
  x = c_xlab,
  y = "Average length intersection",
  title = paste("do Intersection index vs POD:", name),
  color="do category"
) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"))


add_label <- function(plot, label) {
  arrangeGrob(plot, top = textGrob(label, x = unit(0.05, "npc"), y = unit(0.9, "npc"),
                                   just = c("left", "top"),
                                   gp = gpar(fontsize = 16, fontface = "bold")))
}

accum_do <- grid.arrange(
  add_label(plot_jaccard, "A"),
  add_label(prob_intersect, "B"),
  nrow=2)

ggsave(output_jacccard,
       accum_do,
       dpi=300,
       height=6,
       width=6
)


prob_n_do <- ggplot(
  do_df,
  aes(
    x=value
  )
) +
  geom_point(aes(y=n_do_bait, color = "bait DO terms")) +
  geom_line(aes(y=n_do_bait, color = "bait DO terms")) +
  geom_point(aes(y=n_do_prey, color = "prey DO terms")) +
  geom_line(aes(y=n_do_prey, color = "prey DO terms")) +
  facet_wrap(. ~ limit,
             labeller = labeller(
               limit =
                 c("lower_bound_pod" = "Pmin > POD",
                   "upper_bound_pod" = "Pmax < POD")
             ), scales = "free_x") +
labs(
  x = c_xlab,
  y = "Average N DO-annotations",
  title = paste("DO Intersection index vs POD:", name),
  color="do category"
) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"))

ggsave(output_accumulation,
       prob_n_do,
       dpi=300,
       height=3,
       width=6
)