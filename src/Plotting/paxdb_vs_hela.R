library(ggplot2)
library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
df = read.table("intermidiate_data/pca_hela_pax_db.csv", sep="\t", header = T)
df$year <- as.character(df$year)
df[is.na(df$year), "year"] <- ""
df[,"source"] <- "HeLa"
df[df$HeLa == "False","source"] <- "PaxDB"


enough_samples_year  <- rownames(table(df$year)[table(df$year) > 5])
df %>% filter(year %in% enough_samples_year) -> df

# enough_samples_organ <- rownames(table(df$organ)[table(df$organ) > 5])
# df %>% filter(organ %in% enough_samples_organ) -> df

p1 <- ggplot(
  df,
  aes(x=PC1, y=PC2, colour = HeLa)) +
  geom_point() +
  theme_bw() +
  ggtitle("Paxdb vs 0030") +
  theme(legend.position = "none")

p2 <- ggplot(
  df,
  aes(x=PC1, y=PC2, colour = year)) +
  geom_point() +
  theme_bw() + 
  ggtitle("Year produced") +
  theme(legend.position = "none")

add_label <- function(plot, label) {
  arrangeGrob(plot, top = textGrob(label, x = unit(0.05, "npc"), y = unit(0.9, "npc"),
                                   just = c("left", "top"),
                                   gp = gpar(fontsize = 16, fontface = "bold")))
}

grid.arrange(
  add_label(p1, "A"),
  add_label(p2, "B"),
  nrow=1)
