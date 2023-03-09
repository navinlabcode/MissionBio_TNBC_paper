library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ape)
library(phangorn)
library(cowplot)
library(umap)
library(hdf5r)
library(tapestri.cnv)
library(ggplot2)
library(dplyr)
library(magrittr)
library(ggrepel)
library(Rtsne)

# Mission Bio TNBC Paper Figure 3 and 4

# Figure 3A, E, Figure 4A - UMAP
tn_umap <- umap(x, labels=as.factor(x$Cluster), n_neighbors = 15L, a=1, b=1)
iris.labels = as.numeric(as.matrix(x[,"Cluster"]))
plot.iris(tn_umap, iris.labels)

# x = mutation matrix

plot.iris = function(x, labels,
                     main="UMAP",
                     colors=c("red3", "purple4", "#228B22", "blue"),
                     pad=0.1, cex=0.65, pch=19, add=FALSE, legend.suffix="",
                     cex.main=1, cex.legend=1) {
  
  layout = x
  if (is(x, "umap")) {
    layout = x$layout
  } 
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#000000", lwd=1.5)  
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  
  labels.u = unique(labels)
  legend.pos = "topright"
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomright"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  legend(legend.pos, legend=legend.text,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

# Figure 3B, F, Figure 4B - Heatmap

# NGT Mission Bio mutation matrix clustered by heatmap2
x_hm2 <- heatmap.2(as.matrix(NGT_MB_SNV_matrix), trace = "none", na.rm = TRUE, labRow = FALSE, col =colorRampPalette(c("blue","red3"))(100),margins = c(20,5))

# create ordered object of heatmap2 clustering results
x_hm2_ord <- NGT_LC29[,-c(1:2)][rev(x_hm2$rowInd), x_hm2$colInd]

# create csv of ordered heatmap
write.csv(x, "x.csv")

# annotate clusters -> x$Cluster

# = final mutation matrix

# copy number estimations
# calculated from output depth matrix from Mission Bio pipline using standard workflow from Tapestri.cnv R package
# Refer to Tapestri_CNV.Rmd for more in depth code

# visualize heatmap data
col <- c("2" = "black", "1" = "#486B86", "0" = "cream", "NA" = "grey")
Clusters <- as.character(x$Cluster)

# cluster colors
col_list <- list(Clusters = c("1" ="red3", "2" = "purple4", "3" = "#228B22", "4" = "blue"))

# cluster annotation
row_anno <- rowAnnotation(df = as.data.frame(Clusters), col = col_list)

# allele frequency and copy number annotations
# af = (afs()) - allele frequencies per gene
# cn = (cns() - copy numbers per gene

col_af = colorRamp2(c(0,100), c("white", "#228B22"))
col_cn = colorRamp2(c(0,2,4), c("#1565C0", "white", "red3"))
AF_CN = HeatmapAnnotation(AF = anno_barplot(af), CN = anno_barplot(cn), col = list(AF = col_af), col = list(CN = col_cn))

# visualize using complexheatmap
tn_heatmap <- Heatmap(as.matrix(x), 
                        cluster_columns = F, 
                        cluster_rows = F,
                        name = "Genotype",
                        col = col,
                        show_row_names = F, 
                        show_column_names = T,
                        top_annotation = AF_CN)

draw(row_anno + tn_heatmap)

# Figure 3C, G, Figure4C - CADD analysis
# simple bar ggplots of CADD values per gene

# Figure 3D, H, Figure4AD - Neighbor joining trees
tn_dist <- dist(x)
tn_njtree <- nj(tn_dist)
tn_njtree_ggtree <- ggtree(tn_njtree)

# x = mutation matrix

write.table(exprs(tn_njtree_ggtree), "tn_njtree_ggtree_table.csv", sep =",",col.names=T,row.names=F)
# added cluter annotaion from heatmap matrix -> Cluster column

tn_njtree_ggtreee_plot <- tn_njtree_ggtree_table + geom_tippoint(aes(color=Cluster), size=2)
tn_njtree_ggtreee_plot + theme(legend.position="right") + scale_color_manual(values=c("color1", "color2"))

