library(ggplot2)
library(magrittr)
library(tidyverse)
library(arrow)

PLOT_DIR <- "work_folder/plots/POD/plots"
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# Lazy arrow dataset: only the columns used below are scanned, and
# filter/group_by/count are pushed down so we never materialize the
# full 98M-row table in R.
ds <- open_dataset("work_folder/analysis/POD/undirectional/POD_flat.pq") %>%
  select(uniprot_id_bait, uniprot_id_prey, n_tested, n_observed,
         lower_bound_pod, upper_bound_pod)

non_observed <- ds %>% filter(n_observed == 0)

## ---- Plot 1: non-interaction degree distribution ------------------------
degree_df <- bind_rows(
  non_observed %>% count(gene = uniprot_id_bait) %>% collect(),
  non_observed %>% count(gene = uniprot_id_prey) %>% collect()
) %>%
  group_by(gene) %>%
  summarise(degree = sum(n), .groups = "drop")

p_degree <- ggplot(
  degree_df,
  aes(x = degree)
) +
  geom_histogram(fill = "darkorange", bins = 50) +
  scale_x_log10() +
  labs(
    x = expression(log[10]("Non-interaction degree")),
    y = "Number of proteins",
    title = "Full non-interaction degree distribution"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10))

ggsave(file.path(PLOT_DIR, "non_interaction_degree_distribution.png"),
       p_degree, width = 5, height = 3, dpi = 300)

## ---- Plot 2: n_tested density for non-interactions ---------------------
n_tested_counts <- non_observed %>% count(n_tested) %>% collect()
n_tested_total  <- sum(n_tested_counts$n)
n_tested_df <- n_tested_counts %>% mutate(weight = n / n_tested_total)

p_ntested <- ggplot(
  n_tested_df,
  aes(x = n_tested, weight = weight)
) +
  geom_density(fill = "darkorange", alpha = 0.6) +
  labs(
    x = "Number of test",
    y = "Percentage of edges",
    title = "Density of non-interactions vs times tested"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10))

ggsave(file.path(PLOT_DIR, "n_tested_density_non_interaction.png"),
       p_ntested, width = 5, height = 3, dpi = 300)

## ---- Plot 3: lower/upper bound POD scatter + contour, all pairs --------

pod_bounds_df <- ds %>%
  count(n_tested, n_observed, lower_bound_pod, upper_bound_pod, name = "count") %>%
  collect()

FLOOR_LOWER <- 1e-4 # as the lower POD is very low for non-interactions
pod_bounds_df <- pod_bounds_df %>%
  mutate(
    log_lower = log10(pmax(lower_bound_pod, FLOOR_LOWER)),
    log_upper = log10(upper_bound_pod),
    interaction_status = ifelse(
      n_observed == 0, "Not observed", "Observed"
    )
  )

p_pod_bounds <- ggplot(
  pod_bounds_df,
  aes(x = log_lower, y = log_upper, color = interaction_status)
) +
  geom_point(aes(size = count), alpha = 0.6) +
  scale_color_manual(
    values = c("Not observed" = "darkorange", "Observed" = "steelblue")
  ) +
  scale_size_continuous(trans = "log10") +
  labs(
    x = expression(log[10](Q[2.5])),
    y = expression(log[10](Q[97.5])),
    title = "POD credible interval bounds, all protein pairs",
    size = "N protein pairs",
    color = "Interaction"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10))

ggsave(file.path(PLOT_DIR, "pod_bounds_scatter.png"),
       p_pod_bounds, width = 5, height = 3, dpi = 300)

## ---- Plot 4: median non-interaction degree by n_tested -----------------
gene_n_tested_counts <- bind_rows(
  non_observed %>% count(gene = uniprot_id_bait, n_tested) %>% collect(),
  non_observed %>% count(gene = uniprot_id_prey, n_tested) %>% collect()
) %>%
  group_by(gene, n_tested) %>%
  summarise(n = sum(n), .groups = "drop")

max_n_tested <- max(gene_n_tested_counts$n_tested)

median_degree_df <- gene_n_tested_counts %>%
  complete(gene, n_tested = 1:max_n_tested, fill = list(n = 0)) %>%
  arrange(gene, desc(n_tested)) %>%
  group_by(gene) %>%
  mutate(degree = cumsum(n)) %>%
  ungroup() %>%
  filter(degree > 0) %>%
  group_by(n_tested) %>%
  summarise(median_degree = median(degree), .groups = "drop") %>%
  arrange(n_tested)

p_median_degree <- ggplot(
  median_degree_df,
  aes(x = n_tested, y = median_degree)
) +
  geom_line(color = "darkorange") +
  geom_point(color = "darkorange") +
  labs(
    x = "Minimum tests per protein pair",
    y = "Median non-interaction degree",
    title = "Median non-interaction degree by number of tests"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10))

ggsave(file.path(PLOT_DIR, "median_degree_by_n_tested.png"),
       p_median_degree, width = 5, height = 3, dpi = 300)
