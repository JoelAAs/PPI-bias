library(dplyr)
library(arrow)
library(brms)

df <- read_parquet(snakemake@input$detection) %>%
  mutate(
    A_success = as.integer(ms_detection),
    B_success = as.integer(y2h_detection)
  )

# --- 1. Priors: regularize the very rare successes -------------------------
# Intercept centered near logit(base rate of B). If B's rate is ~0.001,
# logit(0.001) ~= -6.9; if ~0.005, ~= -5.3. Adjust to whichever is the response.
priors <- c(
  prior(normal(0, 5),  class = "b"),          # fixed effects (log-odds)
  prior(normal(-6, 3), class = "Intercept"),  # rare-event baseline
  prior(exponential(1), class = "sd")         # node-effect SD
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

# --- 3. Read off the association -------------------------------------------
summary(fit)                                   # check Rhat ~1, no divergences
odds_ratio <- exp(fixef(fit)["A_success", c("Estimate", "Q2.5", "Q97.5")])  # odds ratio + 95% CrI
#odds_ratio

post <- as_draws_df(fit)
#mean(post$b_A_success > 0)

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