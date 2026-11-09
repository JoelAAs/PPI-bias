## power law testig functions copied from github used martaluc
# only change is from a discrete distribution (degree) to
# the interaction rate of protein (interactions/experiments)
# https://github.com/martaluc/powerlaw-ppi-network/blob/main/functions.R 6/11/2026
# https://doi.org/10.7554/eLife.99951
# ->

library(poweRlaw)
library(ggplot2)
library(data.table)


#########################################################################################################################
# function to check the power-law distribution
# input:
# - data = interaction rate vector
# - t = number of threads
# - myseed = seed used to have same random simulated distribution in bootstrapping
# - n_sims = number of simulations in bootstrapping (bootstrap_p function); 100 is the default suggested by poweRlaw developers
# output: list with the bootstrap_p output (p-value, goodness-of-fit, and bootstraps),
#         the empirical CCDF points and the fitted power-law line for plotting
##########################################################################################################################
check_powerLaw <- function(data, t = 1, myseed = 1, n_sims = 100) {
  m_pl <- conpl$new(data)
  est <- estimate_xmin(m_pl)
  m_pl$setXmin(est)
  print(est)
  bs_p <- bootstrap_p(m_pl, threads = t, seed = myseed, no_of_sims = n_sims)
  print(summary(bs_p$bootstraps$xmin))

  xy <- plot(m_pl, draw = FALSE)
  xy$y <- signif(xy$y, digits = 3)

  return(list(
    bs_p = bs_p,
    xy = xy,
    fit = data.table(lines(m_pl, draw = FALSE), model = "PL")
  ))
}

## <-

plot_distributions <- function(fits, dataset) {
  facet_label <- function(ss_limit) paste0("search space > ", ss_limit)
  ss_levels <- facet_label(names(fits))

  # fits[[ss_limit]][[role]] -> one table with search_space_minimum and role columns
  collect <- function(extract) {
    rbindlist(lapply(names(fits), function(ss_limit) {
      rbindlist(lapply(names(fits[[ss_limit]]), function(role) {
        data.table(
          extract(fits[[ss_limit]][[role]]),
          search_space_minimum = factor(
            facet_label(ss_limit),
            levels = ss_levels
          ),
          role = role
        )
      }))
    }))
  }

  points <- collect(function(f) f$xy)
  fit_lines <- collect(function(f) f$fit)
  stats <- collect(function(f) {
    data.table(
      label = sprintf("PL p = %.2f\nLN p = %.2f", f$bs_p_pl$p, f$bs_p_ln$p)
    )
  })

  scaleFUN <- function(x) signif(x, digits = 3)
  ggplot(points, aes(x = x, y = y)) +
    geom_point(size = 0.5) +
    geom_line(data = fit_lines, aes(color = model)) +
    scale_color_manual(values = c(PL = "red", LN = "green"), name = NULL) +
    geom_text(
      data = stats,
      aes(label = label),
      x = -Inf,
      y = -Inf,
      hjust = -0.1,
      vjust = -0.2,
      size = 3
    ) +
    scale_x_continuous(trans = "log2", labels = scaleFUN) +
    scale_y_continuous(trans = "log2", labels = scaleFUN) +
    facet_grid(role ~ search_space_minimum) +
    labs(
      x = "Detectability (detections / tests)",
      y = "Proteins With detectability \u2265 x",
      title = paste("Detectability distribution test of", dataset)
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    )
}


get_interaction_distribution <- function(
  df_test,
  search_space_limits,
  min_tested = 2
) {
  df_test <- df_test[n_tested >= min_tested, ]

  detection_rate <- function(df, col, ss_limit) {
    df_rate <- df[,
      .(p = sum(n_observed) / sum(n_tested), n_pairs = .N),
      by = .(protein = get(col))
    ]
    df_rate <- df_rate[n_pairs > ss_limit]
    df_rate <- df_rate[p > 0]
    return(df_rate)
  }

  search_detectabilities <- list()
  rho_n_df <- data.table()
  for (ss_limit in search_space_limits) {
    bait_detectability_df <- detection_rate(
      df_test,
      "uniprot_id_bait",
      ss_limit
    )
    prey_detectability_df <- detection_rate(
      df_test,
      "uniprot_id_prey",
      ss_limit
    )

    detectability_df <- merge(
      bait_detectability_df,
      prey_detectability_df,
      by = "protein",
      suffixes = c("_bait", "_prey")
    )
    rho <- cor(
      detectability_df$p_bait,
      detectability_df$p_prey,
      method = "spearman"
    )
    n <- nrow(detectability_df)
    rho_n_df <- rbind(
      rho_n_df,
      data.table(ss_limit = ss_limit, n = n, rho = rho)
    )

    detectability_df[, search_space_minimum := ss_limit]

    search_detectabilities[[as.character(ss_limit)]] <- detectability_df
  }
  search_detectabilities[["bp_cor"]] <- rho_n_df
  return(search_detectabilities)
}

check_lognorm <- function(data, t = 1, myseed = 1, n_sims = 100) {
  m_ln <- conlnorm$new(data)
  est <- estimate_xmin(m_ln)
  m_ln$setXmin(est)
  print(est)
  bs_p <- bootstrap_p(m_ln, threads = t, seed = myseed, no_of_sims = n_sims)
  print(summary(bs_p$bootstraps$xmin))

  return(list(
    bs_p = bs_p,
    fit = data.table(lines(m_ln, draw = FALSE), model = "LN")
  ))
}


check_distributions <- function(detectability, t, n_sim = 500) {
  detectability <- detectability[detectability != 0]
  # Filtered out self-interactions are still part of the searchspace.

  pl_fit <- check_powerLaw(detectability, t = t, myseed = 1234, n_sim)
  lognormal_fit <- check_lognorm(detectability, t = t, myseed = 1234, n_sim)

  return(list(
    bs_p_pl = pl_fit$bs_p,
    bs_p_ln = lognormal_fit$bs_p,
    xy = pl_fit$xy,
    fit = rbind(pl_fit$fit, lognormal_fit$fit)
  ))
}

plot_role_correlations <- function(df) {
  g <- ggplot(
    df,
    aes(
      y = rho,
      x = as.factor(ss_limit),
      color = n,
      group = 1
    )
  ) +
    labs(
      y = "Bait-prey detectability correlation",
      x = "Minimum searchspace",
      color = "N passing proteins"
    ) +
    geom_point() +
    geom_line() +
    theme_bw() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    )
  return(g)
}

main <- function() {
  df_pairs <- setDT(arrow::read_parquet(
    snakemake@input$pod_data,
    col_select = c(
      "uniprot_id_bait",
      "uniprot_id_prey",
      "n_observed",
      "n_tested"
    )
  ))
  search_space_limits <- snakemake@params$search_space_minimum

  search_space_distributions <- get_interaction_distribution(
    df_pairs,
    search_space_limits,
  )
  rm(df_pairs)
  detectability_cor <- search_space_distributions[["bp_cor"]]
  search_space_distributions[["bp_cor"]] <- NULL

  fwrite(
    rbindlist(search_space_distributions),
    snakemake@output$probabilities,
    sep = "\t"
  )

  g <- plot_role_correlations(detectability_cor)
  ggsave(
    snakemake@output$plot_bait_prey_cor,
    g,
    dpi = 300,
    height = 4,
    width = 4
  )
  fits <- lapply(search_space_distributions, function(x) {
    list(
      bait = check_distributions(x$p_bait, t = snakemake@threads),
      prey = check_distributions(x$p_prey, t = snakemake@threads)
    )
  })

  g <- plot_distributions(fits, snakemake@wildcards$dataset)
  ggsave(
    snakemake@output$plot_distributions,
    g,
    dpi = 300,
    height = 4 * 2,
    width = 4 * length(fits)
  )
}

main()
