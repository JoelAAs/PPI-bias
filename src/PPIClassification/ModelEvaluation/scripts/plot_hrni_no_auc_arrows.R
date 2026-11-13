library(tidyverse)

# HRNI -> NO AUC per config, one arrow per embedding model.
# Both AUCs are scored on the HRNI test set; "hrni" = trained with HRNI negatives,
# "no" = trained with random (NO) negatives (same pivot as rule hrni_no_delta_auc).
df_auc <- tibble(path = unlist(snakemake@input$metrics)) %>%
    extract(
        path,
        c("source", "classifier", "network_type", "esm_model"),
        "([^/]+)/([^/]+)/all_metrics_([a-z]+)_([A-Za-z0-9]+)\\.csv$",
        remove = FALSE
    ) %>%
    mutate(
        data = map(path, \(f) {
            read_tsv(
                f,
                col_types = cols(model = "c", roc_auc = "d", .default = "-")
            )
        })
    ) %>%
    unnest(data) %>%
    mutate(
        train_type = if_else(str_detect(model, "-random_model_"), "no", "hrni")
    ) %>%
    extract(
        model,
        c("dataset", "neg_limit", "pos_limit", "test_neg"),
        "^([a-zA-Z0-9]+)_[a-zA-Z]+_limit_([0-9.]+)_poslim_([0-9.]+|all)(?:-random)?_model_[A-Za-z0-9]+_metrics_(hrni|no)\\.txt$"
    )
stopifnot(!is.na(df_auc$dataset))

df_auc <- df_auc %>%
    filter(test_neg == "hrni") %>% # always score against the HRNI test set
    pivot_wider(
        id_cols = c(
            source,
            classifier,
            network_type,
            esm_model,
            dataset,
            pos_limit,
            neg_limit
        ),
        names_from = train_type,
        values_from = roc_auc
    ) %>%
    # character sort puts "all" after the numeric pos_limits
    arrange(source, pos_limit, neg_limit) %>%
    mutate(
        source = factor(
            source,
            levels = unlist(snakemake@config$classification_targets$source)
        ),
        config = fct_inorder(paste(pos_limit, neg_limit, sep = " / "))
    )

pd <- position_dodge(width = 0.6)

# position_dodge() doesn't move xend, so xend copies the dodged x after scaling
p <- ggplot(df_auc, aes(x = config, colour = esm_model, group = esm_model)) +
    geom_segment(
        aes(y = hrni, yend = no, xend = stage(config, after_scale = x)),
        position = pd,
        arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
        linewidth = 0.7
    ) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    geom_point(aes(y = hrni, shape = "HRNI"), position = pd, size = 1.8) +
    # invisible, only adds NO to the shape legend (the arrowhead marks NO)
    geom_point(aes(y = no, shape = "NO"), position = pd, alpha = 0) +
    facet_grid(
        network_type ~ source + dataset,
        scales = "free_x",
        space = "free_x"
    ) +
    scale_colour_manual(values = c(ESM2 = "#1b9e77", ProstT5 = "#7570b3")) +
    scale_shape_manual(values = c(HRNI = 16, NO = 17)) +
    guides(shape = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    labs(
        x = "pos_limit / neg_limit",
        y = "ROC AUC (HRNI test set)",
        colour = NULL,
        shape = "Trained on",
        subtitle = "dot = trained on HRNI negatives, arrowhead = trained on NO negatives"
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top"
    )

ggsave(snakemake@output$plot, p, width = 12, height = 6)
