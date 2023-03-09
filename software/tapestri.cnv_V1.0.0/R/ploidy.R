
if (getRversion() >= "2.15.1") utils::globalVariables(c("chr", "name", "cn_median", "chromosome_genes", "cn", "annot.symbol", "cells", "id"))

#' Identify ploidy for samples
#'
#' @param bcfile.norm Normalised barcode distribution
#' @param method Clustering method
#' @param max_cluster number of clusters
#' @param clust barcode clusters based on genotypes
#' @param norm.cluster cluster id for normal cluster
#' @return ploidy
#' @export
#' @examples
#' \dontrun{
#' ploidy <- compute_ploidy(bcfile.norm, "kmeans", max_cluster)
#' ploidy <- compute_ploidy(bcfile.norm, "kmeans", max_cluster, clust, 1)
#' }
#'
compute_ploidy <- function(bcfile.norm, method, max_cluster, clust = NULL, norm.cluster = NULL) {
  barcodes.med <- apply(bcfile.norm, 2, stats::median)
  bcfile.norm <- bcfile.norm[, barcodes.med > 0.0005]

  rowMedian <- function(x, na.rm = FALSE)
    apply(x, 2, stats::median, na.rm = na.rm)

  norm.group.cells <- c()

  if (is.null(clust) | is.null(norm.cluster)) {
    clust <- pick_baseline(bcfile.norm, method, max_cluster)
    if (method == "kmeans") {
      norm.cluster <- which.min(clust$withinss)
      norm.group.cells <- clust$cluster == norm.cluster
    }
    else {
      norm.cluster <- which.max(clust$silinfo$clus.avg.widths)
      norm.group.cells <- clust$cluster == norm.cluster
    }
  } else {
    norm.group.cells <- rownames(clust)[clust$cluster == norm.cluster]
  }
  norm.group <- bcfile.norm[norm.group.cells, ]

  norm.median <- rowMedian(norm.group)

  # ploidy.all = 2 * division(bcfile.norm, norm.median)
  norm.ploidy <- 2 * division(norm.group, norm.median)

  ploidy.all <- norm.ploidy
  for (i in c(1:max_cluster)) {
    if (i != norm.cluster) {
      group.cells <- rownames(clust)[clust$cluster == i]
      if (is.null(group.cells)) {
        group.cells <- clust$cluster == i
      }
      group.df <- bcfile.norm[group.cells, ]
      group.ploidy <- 2 * division(group.df, norm.median)
      ploidy.all <- rbind(ploidy.all, group.ploidy)
    }
  }
  return(ploidy.all)
}


#' Plot line plot of CNV values for each variant
#'
#' @param ploidy ploidy
#' @param amplicons amplicon bed file
#' @param outfile name of output PDF for CNV line plot
#' @param clust clustering output from snv heatmap
#' @param cluster_id normal or control cluster
#' @param min.amplicons filter genes with amplicons < min.amplicons
#' @importFrom magrittr %>%
#' @export
#' @examples
#' \dontrun{
#' ploidy_lineplot(ploidy, amplicons, "lineplot.pdf", clust, 1, 3)
#' }
ploidy_lineplot <- function(ploidy, amplicons, outfile=NULL, clust=NULL, cluster_id = NULL, min.amplicons = 3) {
  if (is.null(ploidy) | is.null(amplicons)) {
    stop("Ploidy or amplicons input is missing")
  }

  cat("Merging barcode and amplicon data \n")
  amplicons <- amplicons[amplicons$name %in% colnames(ploidy), ]
  amplicon_regions <- regioneR::toGRanges(amplicons)
  amplicon_regions <- GenomeInfoDb::sortSeqlevels(amplicon_regions)
  amplicon_regions <- sort(amplicon_regions)

  col_names <- colnames(ploidy)
  filtered.amps <- intersect(amplicon_regions$name, col_names)
  if (length(filtered.amps) == 0) {
    stop("Amplicons do not match with barcode data.")
  }
  col_names <- col_names[!is.na(col_names[match(amplicon_regions$name, col_names)])]

  ploidy.melt <- reshape2::melt(as.matrix(ploidy))
  colnames(ploidy.melt) <- c("cells", "amplicon", "cn")

  if (!is.null(cluster_id) & !is.null(clust)) {
    group.cells <- rownames(clust)[clust$cluster == cluster_id]
    if (is.null(group.cells)) {
        group.cells <- clust$cluster == cluster_id
    }
    group.ploidy <- ploidy[group.cells, ]
    ploidy.melt <- reshape2::melt(as.matrix(group.ploidy))
    colnames(ploidy.melt) <- c("cells", "amplicon", "cn")
  }

  # Select annotations for intersection with regions
  # Note inclusion of custom annotation, and use of shortcuts
  annots <- c("hg19_basicgenes")

  # Build the annotations (a single GRanges object)
  annotations <- annotatr::build_annotations(genome = "hg19", annotations = "hg19_basicgenes")

  dm_annotated <- annotatr::annotate_regions(
    regions = amplicon_regions,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE
  )

  df_dm_annotated <- data.frame(dm_annotated)

  # merge data frames by amplicon ids using "name" column in amplicons file & "amplicon" column in barcodes file)
  ploidy_to_genes <- merge(df_dm_annotated, ploidy.melt, by.x = "name", by.y = "amplicon")
  if (nrow(ploidy_to_genes) == 0 | ncol(ploidy_to_genes) == 0) stop("Cannot merge barcode and amplicon files. Probable reason: panel file does not match barcodes file")

  ploidy_to_genes$name <- as.character(ploidy_to_genes$name)
  ploidy_to_genes <- ploidy_to_genes[!is.na(ploidy_to_genes$annot.symbol), ]
  ploidy_to_genes <- ploidy_to_genes %>%
    dplyr::group_by(cells, annot.symbol) %>%
    dplyr::filter(dplyr::n_distinct(name) >= min.amplicons)

  ploidy_to_genes <- unique(ploidy_to_genes[, c("name", "seqnames", "start", "end", "annot.start", "annot.symbol", "cells", "cn")])
  ploidy_to_genes <- ploidy_to_genes %>%
    dplyr::group_by(annot.symbol) %>%
    dplyr::mutate(cn_median = stats::median(cn, na.rm = T))

  ploidy_to_genes <- as.data.frame(ploidy_to_genes)
  chrom_order <- c(
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
    "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
    "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
    "chr22", "chrX", "chrY", "chrM"
  )

  ploidy_to_genes$seqnames <- factor(x = ploidy_to_genes$seqnames, levels = chrom_order)

  ploidy_to_genes <- ploidy_to_genes[order(ploidy_to_genes$seqnames, ploidy_to_genes$annot.start), ]
  ploidy_to_genes <- ploidy_to_genes %>%
    dplyr::mutate(chromosome_genes = paste(ploidy_to_genes$seqnames, ploidy_to_genes$annot.symbol, sep="_"))
  #ploidy_to_genes.chr_genes <- tidyr::unite_(ploidy_to_genes, "chromosome_genes", c("seqnames", "annot.symbol"))

  ploidy_to_genes$chromosome_genes <- factor(x = ploidy_to_genes$chromosome_genes, levels=unique(ploidy_to_genes$chromosome_genes))

  cat("Creating line plot \n")
  if (!is.null(outfile)) {
    grDevices::pdf(outfile)
  }
  lineplot <- ggplot2::ggplot(data = ploidy_to_genes, ggplot2::aes(y = cn_median, x = chromosome_genes)) +
    ggplot2::stat_summary(
      fun.y = stats::median,
      fun.ymin = stats::median,
      fun.ymax = stats::median,
      geom = "crossbar", width = 0.8
    ) +
    ggplot2::scale_y_continuous(breaks = seq(0, 7, 1), limits = c(0, 7)) +
    ggplot2::theme_bw() + ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      # Change axis line
      axis.line = ggplot2::element_line(colour = "black"),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 5)
    )
  print(lineplot)
  if (!is.null(outfile)) {
    grDevices::dev.off()
  }
  base::rm(ploidy_to_genes)
  devnull <- base::gc()
  return(lineplot)


  # cat("Creating box plot \n")
  # grDevices::pdf(outfile_boxplot)
  # boxplot <- ggplot2::ggplot(data = ploidy_to_genes, ggplot2::aes(y = cn, x = chromosome_genes)) +
  #  ggplot2::geom_boxplot(outlier.shape = NA) +
  #  ggplot2::scale_y_continuous(breaks = seq(0,7, 1), limits = c(0,7)) +
  #  ggplot2::theme(
  #    panel.border = ggplot2::element_blank(),
  #    panel.grid.major = ggplot2::element_blank(),
  #    panel.grid.minor = ggplot2::element_blank(),
  #    # Change axis line
  #    axis.line = ggplot2::element_line(colour = "black"),
  #    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 5)
  #  )
  # print(boxplot)
  # grDevices::dev.off()
}


#' Save heatmap
#'
#' @param ploidy ploidy data
#' @param amplicons amplicons bed records
#' @param outfile outfile to save plots
#' @param groupby.genes group ploidy values by genes
#' @param min.amplicons remove genes with less than min.amplicons number of amplicons
#' @importFrom magrittr %>%
#' @export
#' @examples
#' \dontrun{
#' ploidy_heatmap(ploidy, amplicons, "ploidy.pdf")
#' ploidy_heatmap(ploidy, amplicons, "ploidy.pdf", groupby.genes=TRUE, min.amplicons=3)
#' }
ploidy_heatmap <- function(ploidy, amplicons, outfile=NULL, groupby.genes = FALSE, min.amplicons = 3) {
  if (is.null(ploidy) | is.null(amplicons)) {
    stop("Ploidy or amplicons input is missing")
  }

  amplicon_regions <- regioneR::toGRanges(amplicons)
  amplicon_regions <- GenomeInfoDb::sortSeqlevels(amplicon_regions)
  amplicon_regions <- sort(amplicon_regions)
  col_names <- colnames(ploidy)
  filtered.amps <- intersect(amplicon_regions$name, col_names)
  if (length(filtered.amps) == 0) {
    stop("Amplicons do not match with barcode data.")
  }
  col_names <- col_names[!is.na(col_names[match(amplicon_regions$name, col_names)])]
  ploidy <- ploidy[, filtered.amps]
  chr_names <- GenomeInfoDb::seqnames(amplicon_regions[match(filtered.amps, amplicon_regions$name)])

  column_names_rotation <- 90
  column_names_gp_size <- grid::gpar(fontsize = 4)
  column_split <- chr_names
  if (groupby.genes == TRUE) {
    ploidy.all.df <- reshape2::melt(as.matrix(ploidy))
    colnames(ploidy.all.df) <- c("cells", "amplicon", "cn")
    df_dm_annotated <- annotate_regions_with_genes(amplicon_regions[amplicon_regions$name %in% filtered.amps ], colnames(ploidy))
    # df_dm_annotated$annot.symbol[is.na(df_dm_annotated$annot.symbol)] <- df_dm_annotated$id[is.na(df_dm_annotated$annot.symbol)]
    ploidy.all.df <- merge(df_dm_annotated, ploidy.all.df, by.x = "id", by.y = "amplicon")
    # ploidy.all.df$annot.symbol <- ifelse(is.na(ploidy.all.df$annot.symbol), ploidy.all.df$id, ploidy.all.df$annot.symbol)
    ploidy.all.df <- ploidy.all.df[!is.na(ploidy.all.df$annot.symbol), ]
    ploidy.all.df <- ploidy.all.df %>%
      dplyr::group_by(cells, annot.symbol) %>%
      dplyr::filter(dplyr::n_distinct(id) >= min.amplicons)
    ploidy.all.df <- ploidy.all.df %>%
      dplyr::group_by(cells, annot.symbol) %>%
      dplyr::mutate(cn_median = stats::median(cn, na.rm = T))
    ploidy.all.df <- unique(ploidy.all.df[, c("cells", "annot.symbol", "cn_median")])
    ploidy <- as.data.frame(reshape2::acast(ploidy.all.df, cells ~ annot.symbol, value.var = "cn_median", fun.aggregate = function(x) stats::median(x)))
    sorted_genes <- factor(colnames(ploidy), levels = unique(df_dm_annotated[, 2]))
    ploidy <- ploidy[, as.character(sort(sorted_genes))]
    chr_names <- colnames(ploidy)
    ploidy <- as.matrix(ploidy)
    column_names_rotation <- 0
    column_names_gp_size <- grid::gpar(fontsize = 6, fontface = "bold")
    column_split <- FALSE
    base::rm(ploidy.all.df)
    devnull <- base::gc()
  }
  col_fun <- circlize::colorRamp2(seq(0, 4, length = 3), c("red", "#EEEEEE", "blue"))
  legend_params <- list(title = "Ploidy", col = col_fun, border = "black", title_gp = grid::gpar(fontsize = 8, fontface = "bold"), labels_gp = grid::gpar(fontsize = 6), grid_height = grid::unit(4, "mm"), grid_width = grid::unit(8, "mm"))

  z2 <- scales::squish(ploidy, c(0, 4))
  h <- ComplexHeatmap::Heatmap(data.matrix(z2),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    col = col_fun, name = "cn",
    show_column_dend = FALSE,
    column_names_centered = TRUE,
    heatmap_legend_param = legend_params,
    show_row_dend = FALSE,
    show_row_names = FALSE,
    column_gap = grid::unit(2, "mm"),
    column_title_gp = grid::gpar(fontsize = 10),
    column_names_gp = column_names_gp_size,
    column_names_rot = column_names_rotation,
    use_raster = TRUE,
    raster_quality = 2
  )

  if (!is.null(outfile)) {
    grDevices::pdf(outfile, 15, 7)
  }
  ComplexHeatmap::draw(h)
  if (!is.null(outfile)) {
    grDevices::dev.off()
  }
  return(h)
}

# Divide rows by column median
# @param df barcodes data
# @param column median
# @return df normalized by column median
division <- function(df, z, na.rm = FALSE) {
  df <- as.matrix(df)

  divnorm <- c()

  for (col in colnames(df)) {
    divnorm <- cbind(divnorm, as.numeric(round(df[, col] / z[col], 3)))
  }
  colnames(divnorm) <- colnames(df)
  rownames(divnorm) <- rownames(df)
  return(divnorm)
}
