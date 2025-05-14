library(ggplot2)
library(tidyverse)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)

input   <- args[1]
n_limit <- as.numeric(args[2])
output  <- args[3]

df <- read.table(
  "work_folder/inferred_search_space/analysis/cell_line/bait_wise_prey_plotting.csv",
  sep="\t",
  header=T)

upper_limit = 50
lower_limit = 1/50
p = ggplot(
  df,
  aes(
    x=log10(odds_ratio),
    y=-log10(p_value_adjusted),
    color=cl_id,
    alpha = !(odds_ratio > lower_limit & odds_ratio < upper_limit) 
    )
) +
  geom_point() +
  facet_wrap(cl_id ~.) +
  theme_bw() +
  geom_vline(xintercept = log10(lower_limit), linetype = "dotted") +
  geom_vline(xintercept = log10(upper_limit), linetype = "dotted") +
  theme(legend.position = "none") +
  xlab("log10 OR") +
  ylab("-log10 FDR adjusted p-value")

ratio = 1.618
height = 70
ggsave(output, p, dpi=300, height = height, width =height*ratio, units = "mm")


