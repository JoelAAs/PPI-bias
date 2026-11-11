## power law testig functions copied from github used martaluc
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
  m_pl <- displ$new(data)
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

check_lognorm <- function(data, t = 1, myseed = 1, n_sims = 100) {
  m_ln <- dislnorm$new(data)
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

check_degree_distributions <- function(degree, t, n_sim = 500) {
  degree <- degree[degree != 0]
  # Filtered out self-interactions are still part of the searchspace.

  pl_fit <- check_powerLaw(degree, t = t, myseed = 1234, n_sim)
  # LN needs >= 4 distinct degrees for estimate_xmin; drop the fit if it can't be estimated
  lognormal_fit <- tryCatch(
    check_lognorm(degree, t = t, myseed = 1234, n_sim),
    error = function(e) list(bs_p = list(p = NA), fit = NULL)
  )

  return(list(
    bs_p_pl = pl_fit$bs_p,
    bs_p_ln = lognormal_fit$bs_p,
    xy = pl_fit$xy,
    fit = rbind(pl_fit$fit, lognormal_fit$fit)
  ))
}


get_degree_distro <- function(
  interaction_df,
  detection_limits,
  column_name,
  directed = TRUE
) {
  min_detection <- c(posterior = 0.1, lower_bound_pod = 0.02)
  interaction_df <- interaction_df[
    get(column_name) >= min_detection[[column_name]]
  ]
  degree <- function(df, value_col, role_col, detection_lim) {
    df <- df[get(value_col) >= detection_lim]
    df_deg <- df[,
      .(degree = .N),
      by = .(protein = get(role_col))
    ]
    df_deg <- df_deg[degree > 0]
    return(df_deg)
  }

  detection_limit_degree <- list()
  rho_n_df <- data.table()
  for (detection_limit in detection_limits) {
    bait_degree_df <- degree(
      interaction_df,
      column_name,
      "uniprot_id_bait",
      detection_limit
    )
    prey_degree_df <- degree(
      interaction_df,
      column_name,
      "uniprot_id_prey",
      detection_limit
    )

    degree_df <- merge(
      bait_degree_df,
      prey_degree_df,
      by = "protein",
      suffixes = c("_bait", "_prey"),
      all = TRUE
    )

    if (directed) {
      both_roles <- degree_df[!is.na(degree_bait) & !is.na(degree_prey)]
      rho <- cor(
        both_roles$degree_bait,
        both_roles$degree_prey,
        method = "spearman"
      )
      n <- nrow(both_roles)
      rho_n_df <- rbind(
        rho_n_df,
        data.table(detection_lim = detection_limit, n = n, rho = rho)
      )
    } else {
      degree_df[is.na(degree_bait), degree_bait := 0L]
      degree_df[is.na(degree_prey), degree_prey := 0L]
      degree_df[, degree := degree_bait + degree_prey]
    }

    degree_df[, detection := detection_limit]

    detection_limit_degree[[as.character(detection_limit)]] <- degree_df
  }
  detection_limit_degree[["bp_cor"]] <- rho_n_df
  return(detection_limit_degree)
}


plot_role_degree_correlations <- function(df, directed, lim_column) {
  if (directed) {
    g <- ggplot(
      df,
      aes(
        y = rho,
        x = as.factor(detection_lim),
        color = n,
        group = 1
      )
    ) +
      labs(
        y = "Bait-prey degree correlation",
        x = paste("Minimum", lim_column),
        color = "N proteins in both roles"
      ) +
      geom_point() +
      geom_line() +
      theme_bw() +
      theme(
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10)
      )
  } else {
    g <- ggplot()
  }
  return(g)
}


plot_distributions <- function(fits, title, lim_column) {
  facet_label <- function(limit) paste0(lim_column, " >= ", limit)
  limit_levels <- facet_label(names(fits))

  collect <- function(extract) {
    rbindlist(lapply(names(fits), function(limit) {
      rbindlist(lapply(names(fits[[limit]]), function(role) {
        data.table(
          extract(fits[[limit]][[role]]),
          detection_limit = factor(facet_label(limit), levels = limit_levels),
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
    facet_grid(role ~ detection_limit) +
    labs(
      x = "Degree x",
      y = "Fraction of proteins with degree >= x",
      title = title
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    )
}


check_and_plot <- function(
  infile,
  lim_column,
  detection_limits,
  directed,
  degree_out_file,
  degree_cor_plot,
  dist_plots
) {
  df_pairs <- setDT(arrow::read_parquet(
    infile,
    col_select = c(
      "uniprot_id_bait",
      "uniprot_id_prey",
      lim_column
    )
  ))

  interaction_degree_distributions <- get_degree_distro(
    df_pairs,
    detection_limits,
    lim_column,
    directed
  )

  rm(df_pairs)
  degree_cor <- interaction_degree_distributions[["bp_cor"]]
  interaction_degree_distributions[["bp_cor"]] <- NULL

  fwrite(
    rbindlist(interaction_degree_distributions),
    degree_out_file,
    sep = "\t"
  )

  g <- plot_role_degree_correlations(degree_cor, directed, lim_column)
  ggsave(
    degree_cor_plot,
    g,
    dpi = 300,
    height = 4,
    width = 4
  )

  fits <- lapply(interaction_degree_distributions, function(x) {
    if (directed) {
      list(
        bait = check_degree_distributions(
          x[!is.na(degree_bait), degree_bait],
          t = snakemake@threads
        ),
        prey = check_degree_distributions(
          x[!is.na(degree_prey), degree_prey],
          t = snakemake@threads
        )
      )
    } else {
      list(
        all = check_degree_distributions(x$degree, t = snakemake@threads)
      )
    }
  })

  g <- plot_distributions(
    fits,
    paste(
      snakemake@wildcards$dataset,
      snakemake@wildcards$network_type,
      lim_column
    ),
    lim_column
  )
  ggsave(
    dist_plots,
    g,
    dpi = 300,
    height = 4 * length(fits[[1]]),
    width = 4 * length(fits)
  )
}


main <- function() {
  nt <- snakemake@wildcards$network_type
  if (nt == "directional") {
    directed <- TRUE
  } else if (nt == "undirectional") {
    directed <- FALSE
  } else {
    stop("unknown network_type: ", nt)
  }

  check_and_plot(
    snakemake@input$pod_data,
    "lower_bound_pod",
    snakemake@params$pod_limits,
    directed,
    snakemake@output$pod_degree,
    snakemake@output$plot_bait_prey_cor_pod,
    snakemake@output$plot_distributions_pod
  )

  check_and_plot(
    snakemake@input$fdr_data,
    "posterior",
    snakemake@params$fdr_limits,
    directed,
    snakemake@output$fdr_degree,
    snakemake@output$plot_bait_prey_cor_fdr,
    snakemake@output$plot_distributions_fdr
  )
}

main()
