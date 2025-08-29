library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)

map_entrez <- function(gene_name_df, gene_column="gene_name") {
  entrez_ids <- mapIds(
    org.Hs.eg.db, keys = gene_name_df[, gene_column],
    column = "ENTREZID", keytype = "SYMBOL")
  stacked = stack(entrez_ids)
  stacked_df = as.data.frame(stacked)
  colnames(stacked_df) <- c("entrez_id", gene_column)
  gene_name_df$entrez_id = stacked_df$entrez_id
  return(gene_name_df)
}

dplot_go <- function(df) {
  df$GeneRatio <- sapply(df$GeneRatio, function(x) eval(parse(text=x)))
  df$significant <- df$qvalue < 0.05
  df <- df %>%
    mutate(Description = fct_reorder(Description, GeneRatio))
  p <- ggplot(
    df,
    aes(y=Description,
        x=dataset,
        color=-log10(qvalue),
        size=GeneRatio,
        shape=significant)
    ) + geom_point()
  return(p)
}

dplot_do <- function(df, color="significant") {
  df$significant <- df$qvalue < 0.05
  df$GeneRatio <- sapply(df$GeneRatio, function(x) eval(parse(text=x)))
  df <- df %>%
    mutate(Description = fct_reorder(Description, GeneRatio))
  p <- ggplot(
    df,
    aes(y=Description,
        x=dataset,
        color=!!sym(color),
        size=GeneRatio)
  ) + geom_point()
  return(p)
}


get_enrich_go  <- function(df, n, val_column="norm_degree") {
  df %>% 
    slice_max(!!sym(val_column), n=n) %>%
    pull(entrez_id) -> genes
  
    enrichGO(
      gene          = genes,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP", 
      pAdjustMethod = "BH",
      qvalueCutoff  = 1,
      readable      = TRUE)@result %>%
    return()
}

get_enrich_do  <- function(df, n, val_column="norm_degree") {
  df %>% 
    slice_max(!!sym(val_column), n=n) %>%
    pull(entrez_id) -> genes
  clusterProfiler::enrichDO(genes, qvalueCutoff = 1, pvalueCutoff = 1) %>% 
    return()
}


get_top_id <- function(enrich_df, n, col="ID", value = "qvalue", significant=T) {
  if (significant) {
    enrich_df %>% filter(!!sym(value) < 0.05) -> enrich_df
  }
  enrich_df %>% 
  slice_min(!!sym(value), n=n) %>%
  pull(!!sym(col)) %>% return()
}

get_gene_doid_overlap <- function(enrich_results, top_genes){
  enrich_results@result$geneID_split <- sapply(enrich_results@result$geneID, function(x) str_split(x, "/"))
  i = 1
  n_doids = rep(0, length(enrich_results@gene))
  for (gene in enrich_results@gene) {
    n_doids[i] = sum(sapply(enrich_results$geneID_split, function(x) gene %in% x))
    i = i + 1
  }
  data.frame(
    entrez_id = enrich_results@gene,
    n_doid = n_doids,
    among_top = sapply(enrich_results@gene, function(x) x %in% top_genes)
  ) %>% return()
}
  
get_norm_degree_delta <- function(df1, df2, merged_df, target_col, col_name, missing = T){
  cols <- c("gene_name", target_col)
  delta = merge(df1[, cols], df2[, cols], by="gene_name", all.x=TRUE, suffixes=c("_1", "_2"))
  if (missing){
    delta[is.na(delta[, paste0(target_col, "_2")]), paste0(target_col, "_2")] <- 0
  }
  delta[, col_name] <- delta[, paste0(target_col, "_1")] - delta[, paste0(target_col, "_2")]
  merged_df = merge(merged_df, delta[, c("gene_name", col_name)])
  return(merged_df)
}

hippie_degree="work_folder/degree/full_hippie.csv"
summed_probability="work_folder/degree/flat_summed.csv"
threshold_1="work_folder/degree/flat_min.1.csv"
threshold_2="work_folder/degree/flat_min.2.csv"

df_hp = read.table(hippie_degree, sep ="\t", header=T) %>% 
  map_entrez() %>% filter(!is.na(entrez_id))
df_sp = read.table(summed_probability, sep ="\t", header=T) %>% 
  map_entrez() %>% filter(!is.na(entrez_id))
df_t1 = read.table(threshold_1, sep ="\t", header=T) %>% 
  map_entrez() %>% filter(!is.na(entrez_id))
df_t2 = read.table(threshold_2, sep ="\t", header=T) %>% 
  map_entrez() %>% filter(!is.na(entrez_id))

normalise <- function(x) x/max(x) 

df_hp$norm_degree <- normalise(df_hp$degree)
df_sp$norm_degree <- normalise(df_sp$mean_prey_degree + df_sp$mean_bait_degree)
df_t1$norm_degree <- normalise(df_t1$degree_prey + df_t1$degree_bait)
df_t2$norm_degree <- normalise(df_t2$degree_prey + df_t2$degree_bait)


a = 0.4
degree_dist <- ggplot(
  df_hp,
  aes(x=log(norm_degree))
) + geom_density(aes(fill="HIPPIE"), alpha=a) +
  geom_density(data=df_sp, aes(fill="SP"), alpha=a) +
  geom_density(data=df_t1, aes(fill="T0.1"), alpha=a) +
  geom_density(data=df_t2, aes(fill="T0.2"), alpha=a) +
  theme_bw()+
  labs(
    fill="Dataset",
    x="Log(Normalised Degree)",
    y="Density"
  )

ggsave("work_folder/plots/degree/Distribution.png", degree_dist,
       width= 4,
       dpi=300)

## GO degree
n_genes = 50
df_hp %>% get_enrich_go(n_genes) -> go_hp
go_hp$dataset <- "HIPPIE"
df_sp %>% get_enrich_go(n_genes) -> go_sp
go_sp$dataset <- "SP"
df_t1 %>% get_enrich_go(n_genes) -> go_t1
go_t1$dataset <- "T0.1"
df_t2 %>% get_enrich_go(n_genes) -> go_t2
go_t2$dataset <- "T0.2"

n_id = 10
go_hp %>% get_top_id(n_id) -> go_selected_hp
go_sp %>% get_top_id(n_id) -> go_selected_sp
go_t1 %>% get_top_id(n_id) -> go_selected_t1
go_t2 %>% get_top_id(n_id) -> go_selected_t2


df_go <- bind_rows(go_hp, go_sp, go_t1, go_t2)
df_go %>%  filter(
  ID %in% c(
    go_selected_hp,
    go_selected_sp,
    go_selected_t1,
    go_selected_t2)
  ) -> df_go_top

go_plot <- dplot_go(df_go_top) + theme_bw() + 
  labs(
    x="Dataset",
    y="DO terms",
    shape="q < 0.05"
  ) +
  theme(axis.text.x = element_text(angle = -90))
ggsave("work_folder/plots/degree/GO_enrichment.png",
       height = 5,
       width = 8,
       go_plot, dpi=300)

## DO degree
n_genes = 50
df_hp %>% get_enrich_do(n_genes) -> do_hp
do_hp@result$dataset <- "HIPPIE"
df_sp %>% get_enrich_do(n_genes) -> do_sp
do_sp@result$dataset <- "SP"
df_t1 %>% get_enrich_do(n_genes) -> do_t1
do_t1@result$dataset <- "T0.1"
df_t2 %>% get_enrich_do(n_genes) -> do_t2
do_t2@result$dataset <- "T0.2"

n_id = 5
do_hp@result %>% get_top_id(n_id) -> do_selected_hp
do_sp@result %>% get_top_id(n_id) -> do_selected_sp
do_t1@result %>% get_top_id(n_id) -> do_selected_t1
do_t2@result %>% get_top_id(n_id) -> do_selected_t2


df_do <- bind_rows(
    do_hp@result,
    do_sp@result,
    do_t1@result,
    do_t2@result
    )
df_do %>%  filter(
  ID %in% c(
    do_selected_hp,
    do_selected_sp,
    do_selected_t1,
    do_selected_t2)
  ) -> df_do_top

do_plot <- dplot_do(df_do_top) + theme_bw() + 
  labs(
    x="Dataset",
    y="DO terms",
    color="q < 0.05"
  ) +
  theme(axis.text.x = element_text(angle = -90))
ggsave("work_folder/plots/degree/DO_enrichment.png",
       do_plot,
       height = 3.2,
       dpi=300)


## Degree Difference
delta_degree <- df_hp[, c("gene_name", "norm_degree")]
delta_degree = get_norm_degree_delta(df_hp, df_sp, delta_degree, "norm_degree", "delta_degree_hp_sp")
delta_degree = get_norm_degree_delta(df_hp, df_t1, delta_degree, "norm_degree", "delta_degree_hp_t1")
delta_degree = get_norm_degree_delta(df_hp, df_t2, delta_degree, "norm_degree", "delta_degree_hp_t2")

delta_degree <- delta_degree %>% map_entrez()

n_genes = 50
delta_degree %>% get_enrich_do(n_genes, "delta_degree_hp_sp") -> do_delta_sp
do_delta_sp@result$dataset <- "SP"
delta_degree %>% get_enrich_do(n_genes, "delta_degree_hp_t1") -> do_delta_t1
do_delta_t1@result$dataset <- "T0.1"
delta_degree %>% get_enrich_do(n_genes, "delta_degree_hp_t2") -> do_delta_t2
do_delta_t2@result$dataset <- "T0.2"

n_id = 10
do_delta_sp@result %>% get_top_id(n_id) -> do_delta_selected_sp
do_delta_t1@result %>% get_top_id(n_id) -> do_delta_selected_t1
do_delta_t2@result %>% get_top_id(n_id) -> do_delta_selected_t2

df_do_delta <- bind_rows(
    do_delta_sp@result,
    do_delta_t1@result,
    do_delta_t2@result)

df_do_delta %>% filter(
  ID %in% c(
    do_delta_selected_sp,
    do_delta_selected_t1,
    do_delta_selected_t2)
  ) -> df_do_delta_top

do_delta_plot <- dplot_do(df_do_delta_top, color="qvalue") + theme_bw() + 
  labs(
    x="Dataset",
    y="DO terms"
  ) +
  theme(axis.text.x = element_text(angle = -90))
ggsave("work_folder/plots/degree/DO_delta_enrichment.png",
       do_delta_plot,
       height=3,
       width = 5,
       dpi=300)


do_delta_t1 %>% get_gene_doid_overlap(do_t1@gene) -> n_doid_gene_t1
do_delta_sp %>% get_gene_doid_overlap(do_sp@gene) -> n_doid_gene_sp
do_delta_t2 %>% get_gene_doid_overlap(do_t2@gene) -> n_doid_gene_t2

n_doid_gene_sp$dataset <- "SP"
n_doid_gene_t1$dataset <- "T0.1"
n_doid_gene_t2$dataset <- "T0.2"

df_top_delta <- do.call(
  rbind,
  list(
    n_doid_gene_sp,
    n_doid_gene_t1,
    n_doid_gene_t2)
  )

df_top_delta <- merge(df_top_delta, df_hp, on="entrez_id", all.x=T)

df_top_delta$gene_name <- mapIds(org.Hs.eg.db, keys=df_top_delta$gene_id, column="SYMBOL", keytype="ENTREZID", multiVals="first")
df_top_delta <- df_top_delta %>%
  mutate(gene_name = fct_reorder(gene_name, -norm_degree))


delta_gene_plot <- ggplot(df_top_delta, aes(y = gene_name, x = dataset, fill = n_doid)) +
  geom_tile(color="black") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Gene Name", y = "Dataset", fill = "n_doid") + 
  labs(
    x="Dataset",
    y="Gene name",
    fill="N doids"
  ) 

ggsave("work_folder/plots/degree/genes_top_delta.png", 
       delta_gene_plot,
       height=6,
       dpi = 300)  
