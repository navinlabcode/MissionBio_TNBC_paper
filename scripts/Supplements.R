library(ineq)
library(Gini)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(umap)

# Mission Bio TNBC Papaer Supplemental Figures

# Figure S1

# #A - Gini Index

# amplicons
tn_gini_amplicons <- sapply(x, function(x) ineq(x, type = "Gini"))

# x = depth matrix
# all amplicons avegared acroos each sample

# cells
tn_gini_cells <- apply(x, 1, function(x) ineq(x, type = "Gini"))

# x = depth matrix
# all cells avegared acroos each sample

# B - Coefficient of Variation

# CV
cv = sd(data) / mean(data) * 100

# amplicons
tn_cv_amplicons <- sapply(x, function(x) sd(x) / mean(x) * 100)

# x = depth matrix
# all amplicons avegared acroos each sample

# cells
tn_cv_cells <- apply(x, 1,  function(x) sd(x) / mean(x) * 100)

# x = depth matrix
# all cells avegared acroos each sample

# C
calculated by Mission Bio pipeline

# D
# simple bar ggplot of doublet percentages per sample calcualted in STAR METHODS
# (percentagw of clusters in each heatmap - (doublets/tumor clusters)

# E
calculated by Mission Bio pipeline

# F-I
# standard parameters of seqtk used to downsample fastq files
# metrics calculated by m ission Bio pipeline

# J-K
# standard parameters of seqtk used to downsample fastq files
# mutation matrix is plotted and clusters defined for detection

col <- c("2" = "black", "1" = "#486B86", "0" = "cream", "NA" = "grey")
tn_clusters <- Heatmap(as.matrix(x, 
                         cluster_columns = F, 
                         cluster_rows = F,
                         name = "Clusters",
                         col = col,
                         show_row_names = F, 
                         show_column_names = T)
                         
# x = mutation_matrix

# Figure S2

# A
# same code as Figure S1A, without averages across samples (each mutation value plotted individually)

# B
# simple bar ggplot of depth matrix for each gene per sample

# Figure S3

tn_umap_5samples <- umap(x, labels=as.factor(x$Cluster), n_neighbors = 15L, a=1.2, b=1.2)
iris.labels = as.numeric(as.matrix(x[,"Cluster"]))
plot.iris(tn_umap_5samples, iris.labels)

# x = mutation matrix of all samples combined (tumor cells only)

# Figure S4

# standard Pyclone 2 parameters for eaxch sample
# input allele frequwncy and estimated copy number

# Figure S5

col <- c("2" = "black", "1" = "#486B86", "0" = "cream", "NA" = "grey")
tn_doublets <- Heatmap(as.matrix(, 
                         cluster_columns = F, 
                         cluster_rows = F,
                         name = "Clusters",
                         col = col,
                         show_row_names = F, 
                         show_column_names = T)
                         
 #x = depth mnatrix

(t.test(x, paired = TRUE)) 
# x = data frame of tumor clusters and doublet clusters paired by mutation
# simple bar ggplots of log2(ratio values) of tumor clusters vs doublet clusters

