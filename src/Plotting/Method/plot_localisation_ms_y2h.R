library(ggplot2)
library(ggExtra)
library(tidyverse)
library(magrittr)
library(ggrepel)

df_localisation_count = read.table(
  "data/localisation/gene_to_location.csv",
  sep="\t",
  header=T
) %>% 
  group_by(localisation) %>% 
  summarise(n = n_distinct(gene_name)) %>% 
  ungroup()
  

df_localisation_diff = read.table(
  "work_folder/inferred_search_space/analysis/localisation/same_localisation_method_diff.csv",
  sep="\t",
  header=T
)

df_localisation_ms = read.table(
  "work_folder/inferred_search_space/analysis/localisation/diff_localisation_ms.csv",
  sep="\t",
  header=T
)

df_localisation_y2h = read.table(
  "work_folder/inferred_search_space/analysis/localisation/diff_localisation_y2h.csv",
  sep="\t",
  header=T
)

df_localisation = merge(
  df_localisation_ms,
  df_localisation_y2h,
  by="localisation",
  suffixes = c("_ms", "_y2h")
)

df_localisation = merge(
  df_localisation,
  df_localisation_count,
  by="localisation")

df_localisation_diff = merge(
  df_localisation_diff,
  df_localisation_count,
  by="localisation")

t = 1.2
df_localisation$label <- ifelse(
  abs(df_localisation$OR_ms - df_localisation$OR_y2h) > t,
  df_localisation$localisation,
  "")


## Plot
max_or <- max(c(df_localisation$OR_ms, df_localisation$OR_y2h))
min_or <- min(c(df_localisation$OR_ms, df_localisation$OR_y2h))
p <- ggplot(
  df_localisation,
  aes(
    x=OR_ms,
    y=OR_y2h
    )
  ) + geom_point(aes(size=n)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_text_repel(aes(label = label), size = 2, color="forestgreen") +
  theme_bw() +
  xlim(min_or, max_or) +
  ylim(min_or, max_or) +
  xlab("OR Localisation MS") +
  ylab("OR Localisation Y2H") +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size=8),
    axis.title.y = element_text(size=8)
    #plot.margin = margin(1, 1, 1, 1, "cm")
  )

ggsave("work_folder/plots/localisation_OR_y2h_ms.png", p, height = 8, width = 8, units = "cm", dpi = 300)

