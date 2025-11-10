library(ggplot2)
library(tidyverse)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)
membrane_file = args[1]
output_png = args[2]

membrane_df = read.table(membrane_file, sep="\t", header=TRUE)

g = ggplot(membrane_df, 
       aes(
         y = interaction(membrane_association,target),
         x = size,
         fill = dataset)
       ) +
  geom_col(
    position = position_dodge()
    
    ) +
  geom_errorbar(
    aes(xmin=ci_2.5, xmax=ci_97.5),
    width=0.4,
    colour="black",
    alpha=0.9,
    size=0.3,
    position = position_dodge(width = 0.9)) +
  facet_wrap(~ selection) +
  scale_y_discrete(
    breaks = c("True.bait", "False.bait", "False.prey", "True.prey"),  # values from interaction factor
    labels = c("Membrane\n Bait", "Non Membrane\n Bait", "Non Membrane\n Prey", "Membrane\n Prey")  # custom tick labels
  ) +
  theme_bw() +
  labs(
    fill="Dataset",
    y= "Membrane association - experimental Role",
    x="Average degree") +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -90))

size = 4
ggsave(
  output_png,
  g,
  width = size*1.6,
  height = size,
  dpi = 300
  )

