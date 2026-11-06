## power law testig functions copied from github used martaluc
# https://github.com/martaluc/powerlaw-ppi-network/blob/main/functions.R 6/11/2026
# https://doi.org/10.7554/eLife.99951
# ->

library(poweRlaw)
library(ggplot2)
library(tidyverse)

check_powerLaw <- function(data, plot = F, plot_name ='', t = 1, xlabel, myseed=1, n_sims=100){
  
  m_pl = displ$new(data)
  est = estimate_xmin(m_pl)
  m_pl$setXmin(est)
  #print(est)
  bs_p = bootstrap_p(m_pl,threads = t, seed = myseed, no_of_sims = n_sims)
  #print(summary(bs_p$bootstraps$xmin))
  p <- bs_p$p
  #plot
  if(plot == T){
    xy <- plot(m_pl)
    xy$y <- signif(xy$y, digits = 3)
    l <- lines(m_pl, col = 'green')
    scaleFUN <- function(x) signif(x,digits = 3)
    ggplot(xy, aes(x = x, y = y)) + geom_point() + geom_line(data = l, color = "red") +
      scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2', labels = scaleFUN) + xlab(xlabel) + ylab('Frequency') +
      theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10)) + theme_bw()
    ggsave(paste0(plot_name,'.pdf'),height = 6, width = 6,units = 'cm')
  }
  
  return(bs_p)
}

xmin_estimated <- function(data){
  
  m_pl = displ$new(data)
  est = estimate_xmin(m_pl)
  
  return(est$xmin)
}

## <-

get_protein_degrees <- function(df_pairs) {
    bind_rows(
        df_pairs %>% select(protein = uniprot_id_bait, partner = uniprot_id_prey),
        df_pairs %>% select(protein = uniprot_id_prey, partner = uniprot_id_bait)
    ) %>%
        distinct() %>%
        count(protein, name = "degree") %>%
        return()
}


get_degree_distribution <- function(df_test, search_space_limits) {
    search_space_size <- df_test %>%
        groupby("bait") %>%
        summerise(
            search_space_size =n_distict("prey"),
            )

    search_detectabilities = list()
    for (ss_limit in search_space_limits) {
        detecting_protein <- search_space_size %>%
            filter(search_space_size > ss_limit)

        subset_tests <- df_test %>% 
            filter(bait %in% detecting_protein$bait)

        bait_detectability <- subset_tests %>%
            groupby(bait) %>%
            summerise(p = sum(detections)/count())

        colnames(bait_detectability) <- c("protein", "detectability")
        prey_detectability <- subset_tests %>%
            groupby(prey) %>%
            summerise(p = sum(detections)/count())
        colnames(prey_detectability) <- c("protein", "detectability")
        
        detectability = merge(bait_detectability, prey_detectability, on = protein, how="inner", suffixes = c("_bait", "_prey"))
        detectability["detectability"] =  (detectability["detectability_bait"] +detectability["detectability_prey"])/2

        search_detectabilities[ss_limit] = detectability["detectability"]
    }

    return(search_detectabilities)
 }