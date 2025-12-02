library(ggplot2)
library(tidyverse)
library(reshape2)

df_ms <- read.table(
  "work_folder/per_gene/analysis/interfaces/detection_ms.csv",
  sep="\t",
  header=T
)

df_ms$method = "MS"
df_y2h <- read.table(
  "work_folder/per_gene/analysis/interfaces/detection_y2h.csv",
  sep="\t",
  header=T
)
df_y2h$method = "Y2H"

df_detection = bind_rows(
  df_ms,
  df_y2h
)
df_detection$global_detection_ratio <- df_detection$total_observed/df_detection$total_tests
id_cols <- c(
  "Category",
  "method"
)

df_plot <- melt(df_detection[, c(id_cols, "pair_detection_ratio", "global_detection_ratio")], id.vars = id_cols)
interface_levels <- c("large (strong)", "medium", "small (weak)", "Unknown")
df_plot %>%  
  mutate(Category = factor(Category, interface_levels)) -> df_plot

detection_rate <- ggplot(df_plot,
       aes(
         x=Category,
         y=value,
         color=method,
         group=interaction(method, variable),
         shape=variable)
       ) +
  geom_point() + 
  geom_line() +
  theme_bw() +
  labs(
    x="Interface size",
    y="Detection ratio",
    shape="Detection metric",
    color="Detection Method"
  ) +
  scale_color_manual(
    values = c("darkorange", "blue"),
    labels = c("AP/IP - MS", "Y2H")
  ) +
  scale_shape_manual(
    values = c(8, 15), 
    labels = c("Test average", "Protein pair average")
  ) 

ggsave("work_folder/per_gene/plots/interface_method.png",
       detection_rate,
       dpi=300,
       height=4,
       width=6
)
