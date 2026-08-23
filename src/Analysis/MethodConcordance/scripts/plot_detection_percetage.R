library(dplyr)
library(ggplot2)
library(brms)
library(arrow)
df <- read_parquet(snakemake@input$detection) %>%
    mutate(
        A_success = as.integer(ms_detection),
        B_success = as.integer(y2h_detection)
    )

fit <- readRDS(snakemake@input$fit)

# ---- 1. Raw conditional rates from the data (with Wilson CIs) --------------
raw <- df %>%
    group_by(ms = A_success) %>%
    summarise(k = sum(B_success), n = n(), .groups = "drop") %>%
    rowwise() %>%
    mutate(
        p = k / n,
        lo = binom.test(k, n)$conf.int[1],
        hi = binom.test(k, n)$conf.int[2],
        source = "Naive"
    ) %>%
    ungroup()

# ---- 2. Model-adjusted rates at the average protein (from posterior) -------
post <- as_draws_df(fit)
adj_ms0 <- plogis(post$b_Intercept) # P(Y2H+ | MS-)
adj_ms1 <- plogis(post$b_Intercept + post$b_A_success) # P(Y2H+ | MS+)

adj <- tibble(
    ms     = c(0, 1),
    p      = c(median(adj_ms0), median(adj_ms1)),
    lo     = c(quantile(adj_ms0, .025), quantile(adj_ms1, .025)),
    hi     = c(quantile(adj_ms0, .975), quantile(adj_ms1, .975)),
    source = "Hub-adjusted"
)

plot_df <- bind_rows(raw, adj) %>%
    mutate(ms = factor(ms,
        levels = c(0, 1),
        labels = c("MS negative", "MS positive")
    ))

hub_adjusted_concordance_plot <- ggplot(plot_df, aes(ms, p, colour = source, group = source)) +
    geom_line(position = position_dodge(0.3), linewidth = 0.5, alpha = 0.6) +
    geom_pointrange(aes(ymin = lo, ymax = hi),
        position = position_dodge(0.3), size = 0.7
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
    scale_colour_manual(values = c(
        "Naive" = "blue",
        "Hub-adjusted" = "darkorange"
    )) +
    labs(
        x = "MS Detection",
        y = "P(Y2H detection)",
        colour = NULL,
        title = "Y2H detection among MS pairs",
    ) +
    theme_bw() +
    theme(legend.position = "bottom")


re <- ranef(fit)$mmprotein_aprotein_b[, , "Intercept"] |> as.data.frame()
mixed_effect_dist_plot <- ggplot(re, aes(Estimate)) +
    geom_histogram(bins = 60, fill = "#183552") +
    labs(
        x = "Protein random effect (log-odds)",
        y = "Number of proteins",
        title = "Per-protein random effect distribution"
    ) +
    theme_bw()


ggsave(snakemake@output$mixed_effects_dist,
    mixed_effect_dist_plot,
    dpi = 300, height = 4, width = 4
)
ggsave(snakemake@output$hub_adjusted_concordance,
    hub_adjusted_concordance_plot,
    dpi = 300, height = 4, width = 4
)

