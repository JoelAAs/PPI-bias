library(ggplot2)
library(tidyverse)
library(magrittr)
library(ggrepel)
library(stringr)
library(reshape2)
min_point_y2h = 0.00255416
min_point_ms = 0.002556271
bp_y2h = 0.01762
mf_y2h = 0.2039
bp_ms = 0.03242
mf_ms = 0.1383

specific_point = data.frame(
  limit = rep("lower_bound_pod", 4),
  method = c("y2h", "y2h", "ms", "ms"),
  value = c(min_point_y2h, min_point_y2h, min_point_ms, min_point_ms),
  jaccard = c(bp_y2h, mf_y2h, bp_ms, mf_ms),
  GO = c("ji_bp", "ji_mf", "ji_bp", "ji_mf"),
  desc = c("Y2H","Y2H", "MS", "MS")
)
specific_point$limit <- factor(specific_point$limit, levels = c("upper_bound_pod", "lower_bound_pod"))

greater_go_accumulation_ms <- "work_folder/per_gene/analysis/GO/cumulative/POD_ms_jaccard_greater.csv"
lesser_go_accumulation_ms  <- "work_folder/per_gene/analysis/GO/cumulative/POD_ms_jaccard_lesser.csv"
greater_go_accumulation_y2h <- "work_folder/per_gene/analysis/GO/cumulative/POD_y2h_jaccard_greater.csv"
lesser_go_accumulation_y2h  <- "work_folder/per_gene/analysis/GO/cumulative/POD_y2h_jaccard_lesser.csv"


### Input
df_greater_ms = read.table(
  greater_go_accumulation_ms,
  sep="\t",
  header=T
)
df_greater_ms$method = "ms"
df_lesser_ms = read.table(
  lesser_go_accumulation_ms,
  sep="\t",
  header=T
)
df_lesser_ms$method = "ms"
df_greater_y2h = read.table(
  greater_go_accumulation_y2h,
  sep="\t",
  header=T
)
df_greater_y2h$method = "y2h"
df_lesser_y2h = read.table(
  lesser_go_accumulation_y2h,
  sep="\t",
  header=T
)
df_lesser_y2h$method = "y2h"


go_df = bind_rows(
  df_lesser_ms,
  df_greater_ms,
  df_lesser_y2h,
  df_greater_y2h
)
prob  <- function(x) 1/(1+exp(-x))
logit <- function(x) -log(1/x-1)
go_df[go_df$limit == "upper_bound_pod", "value"] = logit(go_df[go_df$limit == "upper_bound_pod", "value"])
go_df$ji_bp = go_df$sum_ji_bp/go_df$non_na_pairs_ji_bp
go_df$ji_mf = go_df$sum_ji_mf/go_df$non_na_pairs_ji_mf
cols <- c(
  "value",
  "method",
  "limit",
  "ji_bp",
  "ji_mf"
)

go_df_t <- melt(
  go_df[,cols],
  id.vars = c("value", "method", "limit"),
  value.name = "jaccard", 
  variable.name = "GO")


go_df_t$limit <- factor(go_df$limit, levels = c("upper_bound_pod", "lower_bound_pod"))

go_jaccard <- ggplot(
  go_df_t,
  aes(
    x=value,
    y=jaccard,
    color=method,
    shape=GO
  )
) +
  geom_point() +
  geom_line() +
  geom_point() +
  geom_line() +
  
  geom_hline(yintercept=0, linetype="dashed") +
  geom_point(data=specific_point, color="black", size=5, shape=2) +
  geom_text_repel(
     data = specific_point,
     aes(label = desc),
     color = "black",
     size = 4,
     nudge_x = -0.1,
     nudge_y = 0.005
   ) +

  facet_wrap(. ~ limit,
             labeller = as_labeller(
               c("lower_bound_pod" = "P[0.025]",
                 "upper_bound_pod" = "logit(P[0.975])"),
               label_parsed
             ), scales="free_x", strip.position = "bottom",) +
  labs(
    color = "Detection method",
    shape = "GO-category",
    x = "",
    y = "Mean Jaccard index"
  ) +
  theme_bw() +
  scale_color_manual(
    values = c("darkorange", "blue"), 
    labels = c("AP/IP - MS", "Y2H")
  ) +
  scale_shape_manual(
    values = c(8, 15), 
    labels = c("BP", "MF")
  ) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside")

ggsave("manual_figures/GO_jaccard.png",
       go_jaccard,
       dpi=300,
       height=4,
       width=6
)


#### 
go_df_lower = go_df %>% filter(limit == "lower_bound_pod")
lower_count <- ggplot(go_df_lower,
                      aes(
                        x=value,
                        y=log10(non_na_pairs_ji_bp),
                        color=method
                      )) + 
  geom_line(stat="identity") +
  labs(
    x = "",
    y = "Log10 Pairs"
  ) +
  theme_bw() +
  scale_color_manual(
    values = c("darkorange", "blue"), 
    labels = c("AP/IP - MS", "Y2H")
  ) +
  theme(legend.position = "none") +
  scale_y_continuous(position = "right")

ggsave("manual_figures/cumsum_lower.png",
       lower_count,
       dpi=300,
       height=1.5,
       width=3.13
)

colocalisation_df_upper = colocalisation_df %>% filter(limit == "upper_bound_pod")
upper_count <- ggplot(colocalisation_df_upper,
                      aes(
                        x=value,
                        y=log10(non_na_pairs_match_probability),
                        color=method
                      )) + 
  geom_line(stat="identity") +
  labs(
    x = "",
    y = "Log10 Pairs"
  ) +
  theme_bw() +
  scale_color_manual(
    values = c("darkorange", "blue"), 
    labels = c("AP/IP - MS", "Y2H")
  ) +
  theme(legend.position = "none") 

ggsave("manual_figures/cumsum_upper.png",
       upper_count,
       dpi=300,
       height=1.5,
       width=3.13
)
