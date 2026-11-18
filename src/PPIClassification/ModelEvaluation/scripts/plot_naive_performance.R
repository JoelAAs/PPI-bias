library(tidyverse)

# Naive-split performance, box plots of 6 boxes: positive set (facet) x trained on HRNI / NO negatives (x).
# Each box pools seeds x negative source (colour) x embedding (shape); x marks the mean.
#   auc:           ROC AUC on the HRNI test set
#   auc_embedding: same, with embedding on x and training negatives as shape
#   delta_test: ROC AUC on the HRNI test set minus ROC AUC on the NO test set, per model
#   delta_adj:  ROC AUC minus similarity-adjusted ROC AUC, faceted by test set
df_metrics <- read_tsv(
    snakemake@input$all_metrics,
    col_types = cols(
        positive = "c",
        negative = "c",
        seed = "c",
        esm_model = "c",
        roc_auc = "d",
        roc_auc_adj = "d",
        test_set = "c",
        .default = "-"
    )
) %>%
    mutate(
        train_neg = factor(
            if_else(str_ends(negative, "-no"), "NO", "HRNI"),
            levels = c("HRNI", "NO")
        ),
        neg_source = factor(
            str_remove(negative, "-no$"),
            levels = c("fdr", "pod"),
            labels = c("FDR", "POD")
        ),
        positive = factor(
            positive,
            levels = c("hippie", "fdr", "pod"),
            labels = c("HIPPIE", "FDR", "POD")
        ),
        test_set = factor(
            test_set,
            levels = c("exp", "no"),
            labels = c("HRNI test", "NO test")
        )
    )
stopifnot(
    !is.na(df_metrics$positive),
    !is.na(df_metrics$neg_source),
    !is.na(df_metrics$test_set)
)

plot_title <- paste(
    snakemake@wildcards$network_type,
    snakemake@wildcards$classifier
)

# boxes per x (fill), points coloured by negative source with shape for the other grouping
perf_boxes <- function(df, x, shape, y, x_label, shape_label, y_label, hline) {
    ggplot(df, aes(x = {{ x }}, y = {{ y }})) +
        geom_hline(yintercept = hline, colour = "grey60") +
        geom_boxplot(
            aes(fill = {{ x }}),
            outlier.shape = NA,
            width = 0.6,
            alpha = 0.5
        ) +
        geom_point(
            aes(colour = neg_source, shape = {{ shape }}),
            position = position_jitter(width = 0.15, height = 0, seed = 1),
            size = 1.6
        ) +
        stat_summary(
            fun = mean,
            geom = "point",
            shape = 4,
            size = 3,
            stroke = 1.2
        ) +
        scale_fill_manual(
            values = c(
                HRNI = "#d95f02",
                NO = "#7570b3",
                ESM2 = "#66a61e",
                ProstT5 = "#e6ab02"
            ),
            guide = "none"
        ) +
        scale_colour_manual(values = c(FDR = "#1b9e77", POD = "#e7298a")) +
        scale_shape_manual(
            values = c(ESM2 = 16, ProstT5 = 17, HRNI = 16, NO = 17)
        ) +
        labs(
            title = plot_title,
            x = x_label,
            y = y_label,
            colour = "Negative source",
            shape = shape_label
        ) +
        theme_bw() +
        theme(legend.position = "top")
}

p_auc <- df_metrics %>%
    filter(test_set == "HRNI test") %>%
    perf_boxes(
        train_neg,
        esm_model,
        roc_auc,
        "Trained on",
        "Embedding",
        "ROC AUC (HRNI test)",
        0.5
    ) +
    facet_grid(~positive)

p_auc_embedding <- df_metrics %>%
    filter(test_set == "HRNI test") %>%
    perf_boxes(
        esm_model,
        train_neg,
        roc_auc,
        "Embedding",
        "Trained on",
        "ROC AUC (HRNI test)",
        0.5
    ) +
    facet_grid(~positive)

p_delta_test <- df_metrics %>%
    select(positive, neg_source, train_neg, seed, esm_model, test_set, roc_auc) %>%
    pivot_wider(names_from = test_set, values_from = roc_auc) %>%
    mutate(delta = `HRNI test` - `NO test`) %>%
    perf_boxes(
        train_neg,
        esm_model,
        delta,
        "Trained on",
        "Embedding",
        "ΔROC AUC (HRNI test − NO test)",
        0
    ) +
    facet_grid(~positive)

p_delta_adj <- df_metrics %>%
    mutate(delta = roc_auc - roc_auc_adj) %>%
    perf_boxes(
        train_neg,
        esm_model,
        delta,
        "Trained on",
        "Embedding",
        "ROC AUC − adjusted ROC AUC",
        0
    ) +
    facet_grid(test_set ~ positive)

# svglite isn't installed, use the cairo-based grDevices::svg
ggsave(snakemake@output$auc, p_auc, device = svg, width = 7, height = 4.5)
ggsave(snakemake@output$auc_embedding, p_auc_embedding, device = svg, width = 7, height = 4.5)
ggsave(snakemake@output$delta_test, p_delta_test, device = svg, width = 7, height = 4.5)
ggsave(snakemake@output$delta_adj, p_delta_adj, device = svg, width = 7, height = 7)
