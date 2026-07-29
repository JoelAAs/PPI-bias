library(lme4)
library(tidyverse)

df <- read.table(snakemake@input$hri_hrni_file, sep="\t", header=TRUE)

# ---- adjusted point estimate (crossed approximation) ----
fit <- lmer(coexpr ~ X_pos + (1 | uniprot_id_bait) + (1 | uniprot_id_prey),
            data = df)
summary(fit)
confint(fit, "X_posTRUE", method = "Wald")

# ---- vertex bootstrap: crude diff, dependence-correct CI ----
proteins <- unique(c(df$uniprot_id_bait, df$uniprot_id_prey))
n  <- length(proteins)
ia <- match(df$uniprot_id_bait, proteins)   # endpoint index per row
ib <- match(df$uniprot_id_prey, proteins)
val <- df$coexpr
is1 <- df$X_pos == 1
is0 <- !is1

wdiff <- function(w) {
  sum(w * val * is1) / sum(w * is1) -
  sum(w * val * is0) / sum(w * is0)
}

set.seed(0)
B <- 2000
i = 0
boot <- numeric(B)
for (b in seq_len(B)) {
  print(i)
  i = i +1
  cnt <- rmultinom(1, size = n, prob = rep(1/n, n))[, 1]  # resample node counts
  boot[b] <- wdiff(cnt[ia] * cnt[ib])                     # edge weight = product
}

point <- wdiff(rep(1, nrow(df)))
ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("CRUDE: %.4f  95%% CI [%.4f, %.4f]\n", point, ci[1], ci[2]))