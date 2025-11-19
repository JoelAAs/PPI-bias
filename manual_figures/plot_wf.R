library(ggplot2)
library(ggExtra)
library(tidyverse)

setwd("/home/joel/Documents/projects/PPI-bias/manual_figures")
df <- read.table("tested_ms.csv", sep="\t", header=T)
df$group = 1
df[df$n_tested > 2,"group"] = 2

g <- ggplot(df,
       aes(
         x=n_tested,
         y=n_observed,
         color=log10(size)
       )) +
  geom_point(size=2) +
  geom_abline(intercept = 0.5, slope=1)+
  theme_bw() +
  lims(
    x=c(-1,41),
    y=c(-1,41)
  ) +
  labs(
    x="Times tested",
    y="Times interaction observed",
    color="Log10 count pairs"
  ) + theme(
    legend.position = "bottom"
  ) +
  scale_color_continuous(high="darkorange", low="blue")

ggsave("Dot_tested_observed.png",
       g,
       dpi=300,
       height=80,
       width=80,
       units="mm")

g_obs <- 
  df %>% 
  group_by(n_observed) %>% 
  summarise(count=sum(size)) %>% ungroup()

gho = ggplot(
  g_obs,
  aes(
    x=n_observed,
    y=log10(count))) +
  geom_bar(
    stat="identity", fill="darkorange") +
  xlim(-1,41) +
  labs(
    x="",
    y=""
  ) + theme_bw() + coord_flip()

ggsave("hist_observed.png",
       gho,
       dpi=300,
       height=70,
       width=35/1.6,
       units="mm")


g_test <- 
  df %>% 
  group_by(n_tested) %>% 
  summarise(count=sum(size)) %>% ungroup()

ght = ggplot(
  g_test,
  aes(
    x=n_tested,
    y=log10(count))) +
  geom_bar(
    stat="identity", fill="darkorange") +
  xlim(-1,41) +
  labs(
    x="",
    y=""
  ) + theme_bw()

ggsave("hist_tested.png",
       ght,
       dpi=300,
       width=70,
       height=35/1.6,
       units="mm")

df_lower <- read.table("log_lower.csv", sep="\t", header=T)
df_lower$density = df_lower$size/sum(df_lower$size)
df_upper <- read.table("log_upper.csv", sep="\t", header=T)
df_upper$density = df_upper$size/sum(df_upper$size)

prob <- function(x) 1/(1+exp(-x))
df_lower$p <- sapply(df_lower$log_lower, prob)
df_lower$group <- "excluded" 
df_lower$keep <- "keep" 
df_lower[df_lower$p > 0.2, "group"] <- "included" 
df_lower$cum_lower <- sapply(df_lower$p, function(x) sum(df_lower[df_lower$p >= x, "size"]))
df_lower$density <- df_lower$cum_lower/sum(df_lower$cum_lower)

hci_plot <- ggplot(df_lower,
       aes(
         y=log10(cum_lower),
         x=p,
         fill=group
       )) + geom_density(stat="identity", alpha =0.5) +
  geom_vline(
             aes(xintercept=.2), 
             colour="black", linetype="dashed", size=1) +
  scale_fill_manual(values=c("lightblue", "blue")) +
  labs(
    x="2.5 % CI Probability of detection",
    y="Log10 Cumsum protein pairs"
  ) + 
  theme_bw() +
  theme(legend.position = "none")

ggsave("HCI_dist.png",
        hci_plot,
        dpi=300,
        height=80,
        width=80,
        units="mm")
