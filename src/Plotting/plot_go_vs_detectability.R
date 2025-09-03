library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)



### Abundance
ra_pod_go_df = read.table("work_folder/analysis/GO/ra_pod_vs_go_terms.csv", sep="\t", header=T)
ra_pod_go_df_long <- melt(
  ra_pod_go_df, 
  id.vars = c("gene_name", "relative_abundance", "pod"),
  variable.name = "GO_category",
  value.name = "go_count"
  )

ra_pod_go_df_long %>%  
  group_by(GO_category) %>% 
  summarise(
    r_ra = cor(relative_abundance, go_count, method="spearman", use = "complete.obs"),
    r_pod = cor(pod, go_count, method="spearman", use = "complete.obs")
  ) %>% 
  ungroup() -> go_ro

go_ra <- ggplot(
  ra_pod_go_df_long,
  aes(
    x=relative_abundance,
    y=log10(go_count),
    color=GO_category
  )
) + 
  geom_point() + 
  facet_wrap(. ~ GO_category, ncol=1,
             labeller = labeller(
               GO_category =
                 c("n_bp" = "GO: BP",
                   "n_cc" = "GO: CC",
                   "n_mf" = "GO: MF"))) +
  geom_text(data=go_ro, aes(x=2.6, y=2.1, label=paste("ρ:", round(r_ra,3))),
            size=3.5, color="black") +
  theme_bw() +
  labs(
    y="Log10(N GO terms)",
    x="Relative abundance"
  ) +
  theme(legend.position = "none")

go_pod <- ggplot(
  ra_pod_go_df_long,
  aes(
    x=pod,
    y=log10(go_count),
    color=GO_category
  )
) + 
  geom_point() + 
  facet_wrap(. ~ GO_category, ncol=1,
             labeller = labeller(
               GO_category =
                 c("n_bp" = "GO: BP",
                   "n_cc" = "GO: CC",
                   "n_mf" = "GO: MF"))) +
  geom_text(data=go_ro, aes(x=0.35, y=2.1, label=paste("ρ:", round(r_pod,3))),
            size=3.5, color="black") +
  theme_bw() +
  labs(
    y="Log10(N GO terms)",
    x="Rate of detection"
  ) +
  theme(legend.position = "none")


add_label <- function(plot, label) {
  arrangeGrob(plot, top = textGrob(label, x = unit(0.05, "npc"), y = unit(0.9, "npc"),
                                   just = c("left", "top"),
                                   gp = gpar(fontsize = 16, fontface = "bold")))
}



n_go_abundance <- grid.arrange(
  add_label(go_ra, "A"),
  add_label(go_pod, "B"),
  nrow=1)



ggsave("work_folder/plots/GO/ra_pod_n_go.png",
       n_go_abundance,
       dpi=300,
       height=4,
       width=6
       )

#### studies

ns_go_df = read.table("work_folder/analysis/GO/n_studies_go_terms.csv", sep="\t", header=T)
ns_go_df_long <- melt(
  ns_go_df, 
  id.vars = c("gene_name", "n_studies"),
  variable.name = "GO_category",
  value.name = "go_count"
)

ns_go_df_long %>%  
  group_by(GO_category) %>% 
  summarise(
    r_ns = cor(n_studies, go_count, method="spearman", use = "complete.obs"),
  ) %>% 
  ungroup() -> go_ns_cor

go_ns <- ggplot(
  ns_go_df_long,
  aes(
    x=log10(n_studies +1),
    y=log10(go_count),
    color=GO_category
  )
) + 
  geom_point() + 
  facet_wrap(. ~ GO_category, ncol=3,
             labeller = labeller(
               GO_category =
                 c("n_bp" = "GO: BP",
                   "n_cc" = "GO: CC",
                   "n_mf" = "GO: MF"))) +
  geom_text(data=go_ns_cor, aes(x=2.1, y=0.5, label=paste("ρ:", round(r_ns,3))),
            size=3.5, color="black") +
  theme_bw() +
  labs(
    y="Log10(N GO terms + 1)",
    x="Log10(N_studies)"
  ) +
  theme(legend.position = "none")



ggsave("work_folder/plots/GO/n_studies_n_go.png",
       go_ns,
       dpi=300,
       height=2,
       width=6
)

