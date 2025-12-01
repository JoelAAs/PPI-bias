library(ggplot2)
library(tidyverse)
library(GOfuncR)

args <- commandArgs(trailingOnly = TRUE)
  

compare_data = args[1]
output = args[2]

mf_df = read.table(compare_data, sep="\t", header=T)
mf_df$mf_terms <- factor(mf_df$mf_terms, levels = sort(unique(mf_df$mf_terms)))
terms_0.05 = mf_df %>%
  filter(mf_frequency > 0.1)
terms_0.05 <- unique(terms_0.05$mf_terms)

mf_df_filtered <- mf_df %>% filter(mf_terms %in% terms_0.05) 
mf_df_filtered$mf_name <- get_names(mf_df_filtered$mf_terms)$go_name
mf_df_filtered <- mf_df_filtered %>% mutate(
  mf_name = fct_reorder(mf_name, mf_frequency, .fun = max, .desc = TRUE)
)

mf_plot <- ggplot(mf_df_filtered,
       aes(y=mf_frequency,
           x=mf_name,
           fill=method)) +
  geom_bar(position="dodge", stat = "identity") +
  facet_wrap(role~., nrow=3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = -90),
    legend.position = "top"
  ) +
  scale_fill_manual(
    values = c("darkorange", "blue")
  ) +
  labs(
    x="MF GO term",
    y="Frequency among PPIs",
    fill="Dataset"
  )

ggsave(output,
       mf_plot,
       dpi=300,
       height=6,
       width=5)
