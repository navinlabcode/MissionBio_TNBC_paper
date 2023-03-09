# Subgroup clusters based on clustering method
#
# @param bcfile.norm Normalised barcode distribution
# @param method Clustering method
# @param k number of clusters
# @return subgroups
# @examples
# \dontrun{
# sub_group <- tapestri.cnv.ploidy.subgroup(bcfile.norm, "kmeans", 5)
# }
tapestri.cnv.ploidy.subgroup <- function(bcfile.norm, method, k) {
  # k-Means clustering
  if (method == "kmeans") {
    bcfile.clust <- stats::kmeans(bcfile.norm, k)
    sub_group <- bcfile.clust$cluster
  }
  else {
    # hierarchical clustering
    bcfile.clust <- factoextra::hcut(bcfile.norm, k)
    sub_group <- stats::cutree(bcfile.clust, k)
  }
  return(sub_group)
}

# Identify cluster groups for each normalized dataset
#
# @param bcfile.norm Normalised barcode distribution
# @param sub_group sub groups
# @param group group number
# @return cluster group
# @examples
# \dontrun{
# group1 <- tapestri.cnv.ploidy.getClusterGroup(bcfile.norm, sub_group, "1")
# group2 <- tapestri.cnv.ploidy.getClusterGroup(bcfile.norm, sub_group, "2")
# }
tapestri.cnv.ploidy.getClusterGroup <- function(bcfile.norm, sub_group, group = 1) {
  stopifnot(!is.null(group))

  bcfile.group <- subset(sub_group, ((sub_group == group)))

  bcfile.norm.group <- merge(bcfile.group, bcfile.norm, by = "row.names")

  rownames(bcfile.norm.group) <- make.names(bcfile.norm.group[, 1])
  bcfile.norm.group <- bcfile.norm.group[, 6:ncol(bcfile.norm.group)]

  return(bcfile.norm.group)
}


# Save clustering output
#
# @param bcfile.norm Normalised barcode distribution
# @param method Clustering method
# @param k number of clusters
# @param outfile outfile to save plots
# @export
# @examples
# \dontrun{
# tapestri.cnv.ploidy.plot(bcfile.norm, "kmeans", 5, "sample-kmeans.pdf")
# tapestri.cnv.ploidy.plot(bcfile.norm, "hcut", 5, "sample-hcut.pdf")
# }
tapestri.cnv.ploidy.plot <- function(bcfile.norm, method, k, outfile) {
  sub_group <- tapestri.cnv.ploidy.subgroup(bcfile.norm, method, k)
  bcfile.vis <- factoextra::fviz_cluster(list(data = bcfile.norm, cluster = sub_group), stand = FALSE, labelsize = 0)
  ggplot2::ggsave(outfile)
}
