library(ggplot2)
library(tidyverse)
library(magrittr)
library(ggrepel)
library(stringr)
### args
args <- commandArgs(trailingOnly = TRUE)

greater_colocalisation_accumulation_ms <- args[1]
lesser_colocalisation_accumulation_ms  <- args[2]
greater_colocalisation_accumulation_y2h <- args[1]
lesser_colocalisation_accumulation_y2h  <- args[2]

min_point_y2h = 0.00255416
min_point_ms = 0.002556271
  
e_y2h = 12113
o_y2h = 12549
rd_y2h = (o_y2h-e_y2h)/e_y2h

e_ms = 56153
o_ms = 75207
rd_ms = (o_ms-e_ms)/e_ms
  
specific_point = data.frame(
  limit = rep("lower_bound_pod", 2),
  method = c("y2h", "ms"),
  value = c(min_point_y2h, min_point_ms),
  rd = c(rd_y2h, rd_ms),
  desc = c("Y2H", "MS")
)
specific_point$limit <- factor(specific_point$limit, levels = c("upper_bound_pod", "lower_bound_pod"))

greater_colocalisation_accumulation_ms <- "work_folder/per_gene/analysis/localisation/study_match_probability/cumulative/POD_ms_localisation_greater.csv"
lesser_colocalisation_accumulation_ms  <- "work_folder/per_gene/analysis/localisation/study_match_probability/cumulative/POD_ms_localisation_lesser.csv"
greater_colocalisation_accumulation_y2h <- "work_folder/per_gene/analysis/localisation/study_match_probability/cumulative/POD_y2h_localisation_greater.csv"
lesser_colocalisation_accumulation_y2h  <- "work_folder/per_gene/analysis/localisation/study_match_probability/cumulative/POD_y2h_localisation_lesser.csv"


### Input
df_greater_ms = read.table(
  greater_colocalisation_accumulation_ms,
  sep="\t",
  header=T
)
df_greater_ms$method = "ms"
df_lesser_ms = read.table(
  lesser_colocalisation_accumulation_ms,
  sep="\t",
  header=T
)
df_lesser_ms$method = "ms"
df_greater_y2h = read.table(
  greater_colocalisation_accumulation_y2h,
  sep="\t",
  header=T
)
df_greater_y2h$method = "y2h"
df_lesser_y2h = read.table(
  lesser_colocalisation_accumulation_y2h,
  sep="\t",
  header=T
)
df_lesser_y2h$method = "y2h"


colocalisation_df = bind_rows(
  df_lesser_ms,
  df_greater_ms,
  df_lesser_y2h,
  df_greater_y2h
)
prob  <- function(x) 1/(1+exp(-x))
logit <- function(x) -log(1/x-1)
colocalisation_df[colocalisation_df$limit == "upper_bound_pod", "value"] = logit(colocalisation_df[colocalisation_df$limit == "upper_bound_pod", "value"])
c_xlab = "logit(POD)                  Probability of detection"

colocalisation_df$limit <- factor(colocalisation_df$limit, levels = c("upper_bound_pod", "lower_bound_pod"))
colocalisation_df$rd <- (colocalisation_df$sum_localisation_match - colocalisation_df$sum_match_probability)/colocalisation_df$sum_match_probability

prob_localisation <- ggplot(
  colocalisation_df,
  aes(
    x=value,
    y=rd,
    color=method
  )
) +
  geom_point() +
  geom_line() +
  geom_point(data=specific_point, color="black", size=5, shape=2) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_text_repel(
    data = specific_point,
    aes(label = desc),
    color = "black",
    size = 4,
    nudge_x = 0.1,
    nudge_y = 0.05
  ) +
  
  facet_wrap(. ~ limit,
             labeller = as_labeller(
                 c("lower_bound_pod" = "P[0.025]",
                   "upper_bound_pod" = "logit(P[0.975])"),
                 label_parsed
             ), scales="free_x", strip.position = "bottom",) +
  labs(
    color = "Detection method",
    x = "",
    y = "Mean RD (O-E)/E"
  ) +
  theme_bw() +
  scale_color_manual(
    values = c("darkorange", "blue"), 
    labels = c("AP/IP - MS", "Y2H")
  ) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10),
        strip.background = element_blank(),
              strip.placement = "outside")

ggsave("manual_figures/colocal_multi_method.png",
       prob_localisation,
       dpi=300,
       height=4,
       width=6
)


#### 
colocalisation_df_lower = colocalisation_df %>% filter(limit == "lower_bound_pod")
lower_count <- ggplot(colocalisation_df_lower,
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
