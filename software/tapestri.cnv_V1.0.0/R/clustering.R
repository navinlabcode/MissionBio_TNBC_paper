#' Identify max clusters
#'
#' @param bcfile.norm Normalised barcode distribution
#' @param method Clustering method
#' @param outfile Save clusters plot to outfile
#' @return max_cluster
#' @export
#' @examples
#' \dontrun{
#' max_cluster <- cluster_barcodes(bcfile.norm, "kmeans")
#' max_cluster <- cluster_barcodes(bcfile.norm, "kmeans", "sample-cluster.pdf")
#' }
cluster_barcodes <- function(bcfile.norm, method, outfile = NULL) {
  set.seed(42)
  if (method == "kmeans") {
    # Silhouette method
    kmeans.clustering <- factoextra::fviz_nbclust(bcfile.norm, stats::kmeans, method = "silhouette") + ggplot2::labs(subtitle = "k-means")
    if (!is.null(outfile)) {
      ggplot2::ggsave(outfile)
    }
    print(kmeans.clustering)
    n.kmeans.clustering <- kmeans.clustering$data
    max_cluster <- as.numeric(n.kmeans.clustering$clusters[which.max(n.kmeans.clustering$y)])
    return(max_cluster)
  } else {
    # Silhouette method
    hcut.clustering <- factoextra::fviz_nbclust(bcfile.norm, factoextra::hcut, method = "silhouette") + ggplot2::labs(subtitle = "hcut")
    if (!is.null(outfile)) {
      ggplot2::ggsave(outfile)
    }
    print(hcut.clustering)
    n.hcut.clustering <- hcut.clustering$data
    max_cluster <- as.numeric(n.hcut.clustering$clusters[which.max(n.hcut.clustering$y)])
    return(max_cluster)
  }
}


# Identify reference cluster
#
# @param bcfile.norm Normalised barcode distribution
# @param method Clustering method
# @param max_cluster number of clusters
# @param outfile plot clustering output
# @return normalized_groups
# @examples
# \dontrun{
# bcfile.norm <- pick_baseline(bcfile.norm, "kmeans", max_cluster, 'clustering.pdf')
# }
pick_baseline <- function(bcfile.norm, method, max_cluster, outfile = NULL) {
  if (method == "kmeans") {
    # Compute k-means
    normalized_groups <- stats::kmeans(bcfile.norm, centers = max_cluster, nstart = 25)
    # Visualize Clustering Results
    kmeans.vis <- factoextra::fviz_cluster(normalized_groups, data = bcfile.norm, labelsize = 0) + ggplot2::labs(subtitle = "k-means")
    if (!is.null(outfile)) {
      ggplot2::ggsave(outfile)
    }
    # control_cluster_cells <- normalized_groups$cluster == which.min(normalized_groups$withinss)
    return(normalized_groups)
  } else {
    # Hierarchical Clustering
    normalized_groups <- factoextra::hcut(bcfile.norm, centers = max_cluster, nstart = 25)
    hcut.vis <- factoextra::fviz_cluster(normalized_groups, data = bcfile.norm, labelsize = 0) + ggplot2::labs(subtitle = "hcut")
    if (!is.null(outfile)) {
      ggplot2::ggsave(outfile)
    }
    #  define control cluster one with max average Silhouette width
    # control_cluster_cells <- normalized_groups$cluster == which.max(normalized_groups$silinfo$clus.avg.widths)
    return(normalized_groups)
  }
}
