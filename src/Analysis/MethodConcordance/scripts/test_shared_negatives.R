library(dplyr)
library(arrow)
library(brms)

df <- read_parquet(snakemake@input$detection) %>%
  mutate(
    A_success = as.integer(ms_detection),
    B_success = as.integer(y2h_detection)
  )

priors <- c(
  prior(normal(0, 5),  class = "b"),
  prior(normal(-6, 3), class = "Intercept"),
  prior(exponential(1), class = "sd")
)

fit <- brm(
  B_success ~ A_success + (1 | mm(protein_a, protein_b)),
  data    = df,
  family  = bernoulli("logit"),
  prior   = priors,
  chains  = 4, cores = 4,
  iter    = 3000, warmup = 1000,
  control = list(adapt_delta = 0.95),
  seed    = 1
)

summary(fit)
odds_ratio <- exp(fixef(fit)["A_success", c("Estimate", "Q2.5", "Q97.5")])

post <- as_draws_df(fit)

saveRDS(fit, snakemake@output$fit)
writeLines(
  c(
    capture.output(summary(fit)),
    "",
    "--- Odds ratio for A_success (exp of log-odds, with 95% CrI) ---",
    capture.output(odds_ratio)
  ),
  snakemake@output$summary
)