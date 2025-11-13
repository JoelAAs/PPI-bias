library(ggplot2)
library(tidyverse)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)
localisation_file = args[1]
output_png = args[2]

local_df = read.table(localisation_file, sep="\t", header=TRUE)
across_zero <- function(ci_lower, ci_upper) {
  return(ci_lower < 0 & ci_upper > 0)
}

local_df$contain_zero <- across_zero(local_df$ci_2.5, local_df$ci_97.5)
keep_localisation <- unique(local_df[!local_df$contain_zero, "localisation"])
local_df$keep <- local_df$localisation %in% keep_localisation
g = ggplot(local_df[local_df$keep,],
  aes(y=reorder(localisation, delta),
      x=delta,
      fill=role
      )) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(
    aes(xmin=ci_2.5, xmax=ci_97.5, group=role),
    width=0.4,
    colour="black",
    alpha=0.9,
    size=0.3,
    position = position_dodge(width = 0.9)) +
  theme_bw() +
  labs(
    fill="Experimental role",
    x = "PPI localisation fraction: HuRI - Bioplex",
    y="") +
  theme(
    legend.position = "bottom")
  
size = 4
ggsave(
  output_png,
  g,
  width = size*1.6,
  height = size,
  dpi = 300
  )
