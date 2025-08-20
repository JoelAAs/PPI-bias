library(ggplot2)
library(magrittr)
library(tidyverse)
library(gridExtra)
library(grid)

df = read.table("work_folder/inferred_search_space/analysis/bias_reduced_ppis/p_estimated_protein_pairs.csv", sep="\t", header = T)
df <- df %>% filter(gene_name_bait != gene_name_prey) # Needs to be fixed

neg_df <- df %>% filter(n_observed==0)
hci_df <- df %>% filter(n_observed!=0)

p_n <- ggplot(
  neg_df,
  aes(
    x = -log10(p)
  )) +
  geom_histogram(fill="darkorange") +
  scale_y_log10() +
  labs(
    x = "-log10(POD)",
    y = "Number of unique protein pairs"
  ) +
  theme_bw() +
  theme(legend.position = "none") 

p_h <- ggplot(
  hci_df,
  aes(
    x = n_observed,
    y = p,
    fill=as.factor(n_observed)
  )) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.8, linetype="dashed") +
  annotate("label", x = 1.5, y = 0.8, label = paste("p:", 0.8), size=3) +
  labs(
    x = "Times interaction observed",
    y = "Probability of observing interaction"
  ) +
  theme_bw() +
  theme(legend.position = "none") 

add_label <- function(plot, label) {
  arrangeGrob(plot, top = textGrob(label, x = unit(0.05, "npc"), y = unit(0.9, "npc"),
                                   just = c("left", "top"),
                                   gp = gpar(fontsize = 16, fontface = "bold")))
}

grid.arrange(
  add_label(p_h, "A"),
  add_label(p_n, "B"),
  nrow=1)
