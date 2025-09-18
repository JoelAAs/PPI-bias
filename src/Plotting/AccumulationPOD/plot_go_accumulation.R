library(ggplot2)
library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)

min_pairs = 200

## Input
args <- commandArgs(trailingOnly = TRUE)

greater_go_accumulation <- args[1]
lesser_go_accumulation  <- args[2]
name                    <- args[3]
output                  <- args[4]

prob  <- function(x) 1/(1+exp(-x))
logit <- function(x) -log(1/x-1)


go_df = bind_rows(
  greater_go_accumulation,
  lesser_go_accumulation
)

go_df$ji_bp <- go_df$sum_ji_bp/go_df$non_na_pairs_ji_bp
go_df$ji_cc <- go_df$sum_ji_cc/go_df$non_na_pairs_ji_cc
go_df$ji_mf <- go_df$sum_ji_mf/go_df$non_na_pairs_ji_mf

go_df$intersect_bp <- go_df$sum_intersect_bp/go_df$non_na_pairs_intersect_bp
go_df$intersect_cc <- go_df$sum_intersect_cc/go_df$non_na_pairs_intersect_cc
go_df$intersect_mf <- go_df$sum_intersect_mf/go_df$non_na_pairs_intersect_mf


### GO
plot_jaccard <- ggplot(
  go_df,
  aes(
    x=value
  )
) +
  geom_point(aes(y = ji_cc, color = "CC Jaccard")) +
  geom_point(aes(y = ji_bp, color = "BP Jaccard")) +
  geom_point(aes(y = ji_mf, color = "MF Jaccard")) +
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


prob_intersect <- ggplot(
  go_df,
  aes(
    x=value
  )
) +
  geom_point(aes(y = intersect_cc, color = "CC intersect")) +
  geom_point(aes(y = intersect_bp, color = "BP intersect")) +
  geom_point(aes(y = intersect_mf, color = "MF intersect")) +
  facet_wrap(. ~ limit,
             labeller = labeller(
               limit =
                 c("lower_bound_pod" = "Pmin > POD",
                   "upper_bound_pod" = "Pmax < POD")
             ), scales = "free_x") +
labs(
  x = "Probability of detection",
  y = "Average length intersection",
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
  add_label(go_jaccard, "A"),
  add_label(prob_intersect, "B"),
  nrow=2)

ggsave(output,
       abundance_GO,
       dpi=300,
       height=6,
       width=6
)
