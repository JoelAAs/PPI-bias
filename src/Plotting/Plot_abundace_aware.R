library(ggplot2)
library(ggExtra)
library(tidyverse)
library(magrittr)
library(reshape2)
library(ggsignif)
library(grid)
library(gridExtra)

mixture_mean <- function(x, samples, prefix="beta_prediction_", suffix="_mean"){
  mm = 0
  for (c in samples){
    value   = as.numeric(x[paste0(prefix, c, suffix)])
    n_tests = as.numeric(x[paste0("n_tested_", c)])
    if (!is.na(value)) {
      mm = mm + value * n_tests/as.numeric(x["total_tested"])
    }
  }
  return(mm)
}

mixture_var <- function(x, samples){
  mv           = 0
  mix_mean     = as.numeric(x["mixture_mean"])
  total_tested = as.numeric(x["total_tested"])
  for (c in samples){
    mu   = as.numeric(x[paste0("beta_prediction_", c, "_mean")])
    sd   = as.numeric(x[paste0("beta_prediction_", c, "_sd")])
    n_tests = as.numeric(x[paste0("n_tested_", c)])
    if (n_tests != 0) {
      mv = mv + (sd^2 + (mu - mix_mean)^2) * n_tests/total_tested
    }
  }
  return(mv)
}


cnames <- c(
  'gene_name_prey',
  'n_tested_CVCL_0030',
  'n_tested_CVCL_0063',
  'n_tested_CVCL_0291',
  'n_observed_CVCL_0030',
  'n_observed_CVCL_0063',
  'n_observed_CVCL_0291',
  "beta_prediction_CVCL_0030_mean",
  "beta_prediction_CVCL_0030_sd",
  "beta_prediction_CVCL_0063_mean",
  "beta_prediction_CVCL_0063_sd",
  "beta_prediction_CVCL_0291_mean",
  "beta_prediction_CVCL_0291_sd",
  "beta_bait_mean",
  "beta_bait_sd",
  "beta_bait_low_ci",
  "beta_bait_high_ci",
  "n_divergences"
)

samples = c("CVCL_0030", "CVCL_0063", "CVCL_0291")
#df = read.table("intermidiate_data/parameters_abundance_intermid.csv", sep="\t", header = F)
df = read.table("work_folder/analysis/abundance_aware/bait_prey_abundance.csv", sep="\t", header = T)
#df = read.table("intermidiate_data/batch_47_parameters.csv", sep="\t", header = F)
#colnames(df) <- cnames

df <- df %>% filter(n_divergences < 1)
df[,"total_observed"] <- df[,'n_observed_CVCL_0030'] + df[,'n_observed_CVCL_0063'] + df[,'n_observed_CVCL_0291']
df[,"total_tested"]   <- df[,'n_tested_CVCL_0030'] + df[,'n_tested_CVCL_0063'] + df[,'n_tested_CVCL_0291']
df[,"rate_of_detection"] <- df[,"total_observed"]/df[,"total_tested"]

df[, "CL_max"] <- apply(df, 1, function(x) c("CVCL_0030", "CVCL_0063", "CVCL_0291")[
  which.max(x[c('n_tested_CVCL_0030', 'n_tested_CVCL_0063', 'n_tested_CVCL_0291')])])

df[, "mixture_mean"] <- apply(df, 1, function(x) mixture_mean(x, samples))
df[, "mixture_var"]  <- apply(df, 1, function(x) mixture_var(x, samples))

df[,"cof_var"] = df[,"beta_bait_sd"]/df[,"beta_bait_mean"]
df[,"rate_of_detection"] <- df$total_observed/df$total_tested
prob <- function(x) {
  return(1/(1+exp(-x)))
}

df[, "max_bait_probability"] <- prob(df$beta_bait_high_ci)
df[, "min_bait_probability"] <- prob(df$beta_bait_low_ci)

df[, "mixture_RA"] <- 2^df$mixture_mean
df[, "mixture_min_RA"] <- 2^(df$mixture_mean-1.96*df$mixture_var)
df[, "mixture_min_RA"] <- 2^(df$mixture_mean-1.96*df$mixture_var)

df[, "min_prob"] <- df[, "mixture_RA"]*df[, "beta_bait_low_ci"]
df[, "max_prob"] <- df[, "mixture_RA"]*df[, "beta_bait_high_ci"]


df_a <- read.table("data/normalised_log_ra.csv", sep="\t", header=T)
df_a <- melt(df_a, id.vars = "cell_line")
df_a$value <- as.numeric(df_a$value)
df_a %>%
  group_by(cell_line, variable) %>%
  summarise(
    means = mean(value, na.rm=T)
  ) -> df_abundance


df_wide <- df_abundance %>%
  pivot_wider(
    names_from = cell_line,
    values_from = means
  )

df_wide$gene_name_prey <- df_wide$variable
df_wide$variable <- NULL

df <- merge(df, df_wide, by="gene_name_prey")
df[, "mixture_mean_flat"] <- apply(df, 1, function(x) mixture_mean(x, samples,"CVCL_", ""))


### Observed vs Abundance
df$observed <- df$total_observed != 0
counts <- df %>%
  filter(total_tested > 0) %>%
  group_by(total_observed) %>%
  summarize(counts = n())

counts$label <-  sapply(counts$counts, function(x) paste("N:", x))

y_max=max(df$mixture_RA)
y_min=min(df$mixture_RA)

mixmean_obs <- ggplot(
  df,
  aes(
    x = as.factor(total_observed),
    y = mixture_RA,
    fill = as.factor(total_observed)
  )) +
  geom_boxplot() +
  labs(
    y="Weighted harmonised mixture abundance",
    x = "Number of observed interactions") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_signif(
    comparisons = list(c("0", "1")), test=wilcox.test,
    map_signif_level = TRUE, textsize = 4) +
  geom_signif(
    comparisons = list(c("0", "2")), test=wilcox.test,
    map_signif_level = TRUE, textsize = 4, y=y_max+0.5) +
  geom_signif(
    comparisons = list(c("0", "3")), test=wilcox.test,
    map_signif_level = TRUE, textsize = 4, y=y_max+1) +
  geom_hline(yintercept = 1, linetype="dashed") +
  theme_bw() +
    geom_text(
      data = counts,
      aes(
        y = min(df$mixture_RA) - 0.4,  # Adjust x position slightly right of max for visibility
        x = total_observed + 0.7,
        label = label
      ),
      size=3.5,
      inherit.aes = FALSE,
      hjust = 0) +
    theme(legend.position = "none")

ggsave("work_folder/plots/abundance/preyabundance_detection.png",
       mixmean_obs,
       dpi=300,
       height=4,
       width = 6)

# Negatome
p = 0.05
detection_lim =-log(1/p-1)
p_neg <- ggplot(
  df %>% filter(total_observed==0),
  aes(
    x=as.factor(total_tested),
    y = max_prob,
    fill=as.factor(total_tested)
  )) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = detection_lim, linetype="dashed") +
  annotate("label", x = 1, y = detection_lim, label = paste("p:", p)) +
  labs(
    x = "Number of tests",
    y = "Upper bound logit(P)"
       ) +
  theme_bw() +
  theme(legend.position = "none") 


### HCI
p = 0.90
detection_lim =-log(1/p-1)
p_hci  <- ggplot(
  df %>% filter(total_observed!=0),
  aes(
    x=as.factor(total_observed),
    y = prob(min_prob),
    fill=as.factor(total_observed)
  )) +
  geom_boxplot() +
  #ylim(-8, 1) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = p, linetype="dashed") +
  annotate("label", x = 1, y = p, label = paste("p:", p)) +
  labs(
    x = "Number of observed",
    y = "Lower bound P"
  ) +
  theme_bw() +
  theme(legend.position = "none") 

add_label <- function(plot, label) {
  arrangeGrob(plot, top = textGrob(label, x = unit(0.05, "npc"), y = unit(0.9, "npc"),
                                   just = c("left", "top"),
                                   gp = gpar(fontsize = 16, fontface = "bold")))
}

hci_neg <- grid.arrange(
  add_label(p_hci, "A"),
  add_label(p_neg, "B"),
  nrow=1)

ggsave("work_folder/plots/abundance/pos_neg.png",
       hci_neg,
       dpi=300,
       height=4,
       width = 8)



### Upper vs Lower abundance modification
ggplot(
  df,
  aes(
    x=min_prob,
    y=max_prob,
    colour = total_tested
  )) +
  geom_point() +
  scale_color_continuous(high="orange" ,low="darkgreen") +
  theme_bw()


### Upper vs Lower abundance modification
p_dist_pos <- ggplot(
  df %>% filter(total_observed != 0),
  aes(
    x=prob(lower_bound_pod),
    y=mixture_mean,
    colour = rate_of_detection
  )) +
  geom_point() +
  scale_color_continuous(high="orange" ,low="darkgreen") +
  theme_bw() +
  theme(legend.position = "top") +
  labs(
    x="Lower bound POD",
    y="Weighted abundance"
  )


p_dist_neg <- ggplot(
  df %>% filter(total_observed == 0),
  aes(
    x=upper_bound_pod,
    y=mixture_mean,
    colour = total_tested
  )) +
  geom_point() +
  scale_color_continuous(high="orange" ,low="darkgreen") +
  theme_bw() +
  theme(legend.position = "top") +
  labs(
    x="Upper bound logit(POD)",
    y="Weighted abundance"
    )

max_min <- grid.arrange(
  add_label(p_dist_neg, "A"),
  add_label(p_dist_pos, "B"),
  nrow=1)


ggsave("work_folder/plots/abundance/max_min.png",
       max_min,
       dpi=300,
       height=4,
       width = 8)
