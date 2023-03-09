library(flowCore)
library(tidyverse)
library(cowplot)
library(copykit)
library(BiocParallel)

# Mission Bio TNBC Paper Figure 2

# A - Flow data
tn_fcs <- read.FCS(x, transformation=FALSE)
write.table(exprs(tn_fcs), "tn_fcs_table.csv", sep =",",col.names=T,row.names=F)
ggplot(tn_fcs_table, aes(x=DAPI-A)) + geom_density(fill="#4472C4")

# x = flow data (fcs file per sample)

# B - Pseudo-bulk data
# standard workflow from Minussi et al. (2021)
# https://navinlabcode.github.io/CopyKit-UserGuide/

# Parallelization
register(MulticoreParam(progressbar = T, workers = 8), default = T)

# Run pre-processing module
tumor <- runVarbin("/path/to/marked/bam/files/",
                 remove_Y = TRUE)

# Mark euploid cells if they exist
tumor <- findAneuploidCells(tumor)

# Mark low-quality cells for filtering
tumor <- findOutliers(tumor)

# Visualize cells labeled by filter and aneuploid status
plotHeatmap(tumor, label = c('outlier', 'is_aneuploid'), row_split = 'outlier')

# Remove cells marked as low-quality and/or aneuploid from the copykit object
tumor <- tumor[,SummarizedExperiment::colData(tumor)$outlier == FALSE]
tumor <- tumor[,SummarizedExperiment::colData(tumor)$is_aneuploid == TRUE]

# kNN smooth profiles
tumor <- knnSmooth(tumor)

# Create a umap embedding 
tumor <- runUmap(tumor)

# Calculate Integer values
tumor <- calcInteger(tumor, method = 'scquantum', assay = 'smoothed_bincounts')

# Search for the K value that maximizes jaccard similarity for clustering of subclones 
# Plot the results
# This step is optional. A fixed K value can be provided later to findClusters()
tumor <- findSuggestedK(tumor)
plotSuggestedK(tumor)

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
tumor  <- findClusters(tumor)
plotUmap(tumor, label = 'subclones')

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
tumor <- calcConsensus(tumor)
tumor <- runConsensusPhylo(tumor)

# Plot a copy number heatmap with clustering annotation
plotHeatmap(tumor, label = 'subclones')

# We use the plotHeatmap function to plot the
# consensus matrix by passing the argument `consensus`.
plotHeatmap(tumor,
            consensus = TRUE,
            label = 'subclones',
            group = 'spatial_info')

# Plot ration and integer values
plotRatio(tumor)
plotRatio(tumor, "cell")

# Plot consensus sequencews
plotConsensusLine(tumor)

# C - Mutation types and signatures
# simple stacked bar ggplots of mutation types and mutation signatures

# D - Bulk exome allele frequencies
# simple bar ggplots of bulk allele frequencies per gene



