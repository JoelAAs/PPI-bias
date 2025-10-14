if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="https://cran.mirror.garr.it/CRAN/")

BiocManager::install("HDO.db")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

install.packages(
    c("tidyverse", "ggplot2", "reshape2"),
     repos="https://cran.mirror.garr.it/CRAN/")
