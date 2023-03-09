library(tidyverse)
library(cowplot)

# Mission Bio TNBC Paper Figure 5

#A - Consensus pseudo-bulk vs Single Cell Mission Bio copy number estimations
# Pearson correlation test
cor.test(x, y,  method = c("pearson"))

# x = Mission Bio copy number estimation
# y = Consensus Pseudo-bulk copy number

# simple dot ggplots of copy numbers (Mission Bio vs Pseudo-bulk)

# B - Bulk exome vs Single Cell Mission Bio allele frequencies
# Pearson correlation test
cor.test(x, y,  method = c("pearson"))
# x = Mission Bio allele frequency
# y = Bulk exome allele frequency

# simple dot ggplots of allele frequencies (Mission Bio vs Bulk)

# C - Pyclone2 vs Mission Bio cluster comparisons
# simple stacked bar ggplots of cluster distributions from Figure S4 (Pyclone2) vs Single Cell Heatmaps (Figure 3B, 3F, 4B)  
