library(tidyverse)

# ROC AUC per config, four lines per config (embedding = colour, test set = linetype),
# each from the AUC trained with HRNI negatives to the AUC trained with random (NO) negatives (shape).
# Drawn twice: on the similarity-adjusted and on the unadjusted AUC.
df_metrics <- read_tsv(
    snakemake@input$metrics,
    col_types = cols(
        source = "c",
        model = "c",
        roc_auc = "d",
        roc_auc_no = "d",
        roc_auc_adj = "d",
        roc_auc_adj_no = "d",
        .default = "-"
    )
) %>%
    mutate(
        train_type = if_else(str_detect(model, "-random_model_"), "no", "hrni")
    ) %>%
    extract(
        model,
        c("dataset", "network_type", "neg_limit", "pos_limit", "esm_model"),
        "^([a-zA-Z0-9]+)_([a-z]+)_limit_([0-9.]+)_poslim_([0-9.]+|all)(?:-random)?_model_([A-Za-z0-9]+)$"
    )
stopifnot(!is.na(df_metrics$dataset))

plot_auc_lines <- function(df_metrics, auc_hrni_test, auc_no_test, y_label) {
    df_auc <- df_metrics %>%
        select(
            source,
            network_type,
            esm_model,
            dataset,
            pos_limit,
            neg_limit,
            train_type,
            "HRNI test" = all_of(auc_hrni_test),
            "NO test" = all_of(auc_no_test)
        ) %>%
        pivot_longer(
            c("HRNI test", "NO test"),
            names_to = "test_set",
            values_to = "auc"
        ) %>%
        mutate(test_set = factor(test_set, levels = c("HRNI test", "NO test"))) %>%
        pivot_wider(names_from = train_type, values_from = auc) %>%
        # character sort puts "all" after the numeric pos_limits
        arrange(source, pos_limit, neg_limit) %>%
        mutate(
            source = factor(
                source,
                levels = unlist(snakemake@config$classification_targets$source)
            ),
            config = fct_inorder(paste(pos_limit, neg_limit, sep = " / "))
        )

    # the trained-on-HRNI and trained-on-NO AUCs as points, for the shape aesthetic
    df_points <- df_auc %>%
        pivot_longer(c(hrni, no), names_to = "train_neg", values_to = "auc") %>%
        mutate(
            train_neg = factor(
                train_neg,
                levels = c("hrni", "no"),
                labels = c("HRNI", "NO")
            )
        )

    # segments and points share the group so they dodge to the same x;
    # position_dodge() doesn't move xend, so xend copies the dodged x after scaling
    pd <- position_dodge(width = 0.7)

    ggplot(df_auc, aes(x = config, colour = esm_model)) +
        geom_hline(yintercept = 0.5, colour = "grey60") +
        geom_segment(
            aes(
                y = hrni,
                yend = no,
                xend = stage(config, after_scale = x),
                linetype = test_set,
                group = interaction(esm_model, test_set)
            ),
            position = pd,
            linewidth = 0.6
        ) +
        geom_point(
            data = df_points,
            aes(
                y = auc,
                shape = train_neg,
                group = interaction(esm_model, test_set)
            ),
            position = pd,
            size = 1.8
        ) +
        facet_grid(
            network_type ~ source + dataset,
            scales = "free_x",
            space = "free_x"
        ) +
        scale_colour_manual(values = c(ESM2 = "#1b9e77", ProstT5 = "#7570b3")) +
        scale_linetype_manual(
            values = c("HRNI test" = "solid", "NO test" = "dashed")
        ) +
        scale_shape_manual(values = c(HRNI = 16, NO = 2)) +
        labs(
            x = "pos_limit / neg_limit",
            y = y_label,
            colour = "Embedding",
            linetype = "Test set",
            shape = "Trained on"
        ) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "top",
            legend.key.width = unit(1.2, "cm")
        )
}

p_adj <- plot_auc_lines(
    df_metrics,
    "roc_auc_adj",
    "roc_auc_adj_no",
    "Similarity-adjusted ROC AUC"
)
p_unadj <- plot_auc_lines(df_metrics, "roc_auc", "roc_auc_no", "ROC AUC")

ggsave(snakemake@output$plot, p_adj, width = 14, height = 6)
ggsave(snakemake@output$plot_unadjusted, p_unadj, width = 14, height = 6)
