library(ggplot2)
library(ggExtra)

df = read.table("tests_models.csv", sep="\t", header = T)

df[,"total_observed"] <- df[,'n_observed_CVCL_0030'] + df[,'n_observed_CVCL_0063'] + df[,'n_observed_CVCL_0291']
df[,"total_tested"] <- df[,'n_tested_CVCL_0030'] + df[,'n_tested_CVCL_0063'] + df[,'n_tested_CVCL_0291']
df[,"rate_of_detection"] <- df[,"total_observed"]/df[,"total_tested"]

df[, "CL_max"] <- apply(df, 1, function(x) c("CVCL_0030", "CVCL_0063", "CVCL_0291")[
  which.max(x[c('n_tested_CVCL_0030', 'n_tested_CVCL_0063', 'n_tested_CVCL_0291')])])

df[,"mixure_mean"] <- df["beta_prediction_0030_mean"]*df["n_tested_CVCL_0030"]/df[,"total_tested"] +
  df["beta_prediction_0063_mean"]*df["n_tested_CVCL_0063"]/df[,"total_tested"] +
  df["beta_prediction_0291_mean"]*df["n_tested_CVCL_0291"]/df[,"total_tested"] 

df[,"mixure_var"] = 
  df["n_tested_CVCL_0030"]/df[,"total_tested"]*(df["beta_prediction_0030_mean"]^2 + df["beta_prediction_0030_sd"]^2) +
  df["n_tested_CVCL_0063"]/df[,"total_tested"]*(df["beta_prediction_0063_mean"]^2 + df["beta_prediction_0063_sd"]^2) +
  df["n_tested_CVCL_0291"]/df[,"total_tested"]*(df["beta_prediction_0291_mean"]^2 + df["beta_prediction_0291_sd"]^2) -
  df[,"mixure_mean"]^2

df[,"cof_var"] = df[,"beta_bait_sd"]/df[,"beta_bait_mean"]
df[df[,"total_observed"] == 0, "total_observed"] = NA
prob <- function(x) {
  return(1/(1+exp(-x)))
}

p <- ggplot(
  df,
  aes(
    y = mixure_mean - 1.96*sqrt(mixure_var),
    x = beta_bait_mean + 1.96*beta_bait_sd,
    colour = total_tested
  )
) + 
  geom_point() +
  geom_density_2d() +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_gradient2(mid = "#12875C", high = "#FF8615") +
  ylab("Lower bounds 95% log(OR detection) untargeted") +
  xlab("Upper bounds 95% log(OR detection) baited")

ggMarginal(p, type = "histogram")

df[, "max_bait_probability"] <- prob(df$beta_bait_mean + 1.96*df$beta_bait_sd)
df[, "min_bait_probability"] <- prob(df$beta_bait_mean - 1.96*df$beta_bait_sd)

ggplot(
  df,
  aes(
    x=log(max_bait_probability),
    y=log(min_bait_probability),
    color=rate_of_detection
  )
) + geom_point() +
  theme_bw() +
  xlab("Upper bounds 95% log(OR detection) baited") +
  ylab("lower bounds 95% log(OR detection) baited") 
  
ggplot(
  df,
  aes(
    y = beta_bait_sd,
    x = mixure_mean,
    colour = CL_max
  )
) + 
  geom_point() +
  theme_bw() 



 ggplot(
  df,
  aes(
    x=rate_of_detection,
    y=log(min_bait_probability),
    color=total_observed
  )
) + geom_point() +
  theme_bw() 
  
