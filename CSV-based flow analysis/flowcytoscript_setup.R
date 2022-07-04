# Installation of required packages for csv version of flowcytoscript
# Run this once to set things up.
# Before running, you may wish to update R and/or RStudio.

install.packages("digest")
install.packages("dunn.test")
install.packages("ggplot2")
install.packages("ggridges")
install.packages("RANN")
install.packages("RColorBrewer")
install.packages("reshape2")
install.packages("Rtsne")
install.packages("umap")
install.packages("dplyr")
install.packages("FNN")
install.packages("devtools")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("FlowSOM")

library(devtools)
devtools::install_github('exaexa/EmbedSOM')