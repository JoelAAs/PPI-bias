library(ggplot2)
library(ggExtra)
library(tidyverse)
library(magrittr)
library(reshape2)
library(ggsignif)


mixture_mean <- function(x, samples, prefix="beta_prediction_", suffix="_mean"){
  mm = 0
  for (c in samples){
    value   = as.numeric(x[paste0(prefix, c, suffix)])
    n_tests = as.numeric(x[paste0("n_tested_CVCL_", c)])
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
    n_tests = as.numeric(x[paste0("n_tested_CVCL_", c)])
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
  "beta_prediction_0030_mean",
  "beta_prediction_0030_sd",
  "beta_prediction_0063_mean",
  "beta_prediction_0063_sd",
  "beta_prediction_0291_mean",
  "beta_prediction_0291_sd",
  "beta_bait_mean",
  "beta_bait_sd",
  "lower_bound_bait",
  "upper_bound_bait",
  "n_divergences"
)

samples = c("0030", "0063", "0291")
df = read.table("intermidiate_data/parameters_abundance_intermid.csv", sep="\t", header = F)
#df = read.table("intermidiate_data/batch_47_parameters.csv", sep="\t", header = F)
colnames(df) <- cnames

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

df[, "max_bait_probability"] <- prob(df$upper_bound_bait)
df[, "min_bait_probability"] <- prob(df$lower_bound_bait)

df[, "min_mixture_RA"]       <- 2^(df$mixture_mean - 1.96*sqrt(df$mixture_var))
df[, "mixture_RA"]           <- 2^(df$mixture_mean)
df[, "max_mixture_RA"]       <- 2^(df$mixture_mean + 1.96*sqrt(df$mixture_var))
df[, "logit_mean"]           <- df$mixture_RA*df$beta_bait_mean

df[, "min_prob"] <- df[, "mixture_RA"]*df[, "lower_bound_bait"]
df[, "max_prob"] <- df[, "mixture_RA"]*df[, "upper_bound_bait"]


df_a <- read.table("data/normalised_log_ra.csv", sep="\t", header=T)
df_a <- melt(df, id.vars = "cell_line")
df$value <- as.numeric(df_a$value)
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

y_max=max(df$mixture_mean)
y_min=mmin(df$mixture_mean)

ggplot(
  df,
  aes(
    x = as.factor(total_observed),
    y = mixture_mean,
    fill = as.factor(total_observed)
  )) +
  geom_boxplot() +
  labs(
    x="Mixture log2 relative abundance",
    y = "Number of observed interactions") +
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

  theme_bw() +
    geom_text(
      data = counts,
      aes(
        y = min(df$mixture_mean) - 0.5,  # Adjust x position slightly right of max for visibility
        x = total_observed + 0.9,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = 0) +
    theme(legend.position = "none")
  

# Negatome
p = 0.05
detection_lim =-log(1/p-1)
ggplot(
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
    y = "Upper CI beta bait at expected abundance"
       ) +
  theme_bw() +
  theme(legend.position = "none") 


### HCI
p = 0.6
detection_lim =-log(1/p-1)
ggplot(
  df %>% filter(total_observed!=0),
  aes(
    x=as.factor(total_observed),
    y = prob(min_prob),
    fill=as.factor(total_observed)
  )) +
  geom_boxplot() +
  ylim(-8, 1) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = detection_lim, linetype="dashed") +
  annotate("label", x = 1, y = detection_lim, label = paste("p:", p)) +
  labs(
    x = "Number of observed",
    y = "Lower CI beta bait at expected abundance"
  ) +
  theme_bw() +
  theme(legend.position = "none") 


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

