# DEPRECATED
# LOH analysis
#
# @param lfile loom file
# @param bc.df normalized barcode distribution data
# @param z.filter enable zygosity filter
# @param z.variants zygosity percentage cutoff
# @param z.missing missing zygosity cutoff
# @param z.het.min minimum heterozygous cutoff
# @param gt.filter enable genotype quality filter
# @param gt.gqc Genotype quality cutoff (default 30)
# @param gt.dpc Read depth cutoff (default 10)
# @param gt.afc Allele frequency cutoff (default 20)
# @param gt.mv Remove variants genotyped in < X percent of cells (default 50)
# @param gt.mc Remove cells with mutations in < X percent of variants (default 50)
# @param gt.mm Remove variants mutated in < X percent of cells (default 1)
# @param gt.mask mask low quality GT as missing (if GQ/DP/AF lower than cutoff, default FALSE)
# @return loh data
# @examples
# \dontrun{
# loh_data <- tapestri.cnv.loh(lfile, bc.df)
# loh_data <- tapestri.cnv.loh(lfile, bc.df,
#   z.filter = TRUE, z.variants = 90, z.missing = 100,
#   z.het.min = 5, gt.filter = FALSE, gt.gqc = 30, gt.dpc = 10, gt.afc = 20, gt.mv = 50, gt.mc = 50, gt.mm = 1,
#   gt.mask = FALSE
# )
# }
tapestri.cnv.loh <- function(lfile, bc.df, z.filter = TRUE, z.variants = 90, z.missing = 100, z.het.min = 5, gt.filter = FALSE, gt.gqc = 30, gt.dpc = 10, gt.afc = 20, gt.mv = 50, gt.mc = 50, gt.mm = 1, gt.mask = FALSE) {

  # vectors with information from rows
  panelannots <- tapestri.cnv.loh.annotations(lfile)

  # matrix with genotypic information
  if (gt.filter == TRUE) {
    f <- tapestri.cnv.loh.genotypes.filter(lfile, gqc = gt.gqc, dpc = gt.dpc, afc = gt.afc, mv = gt.mv, mc = gt.mc, mm = gt.mm, gt.mask = gt.mask)
  } else {
    f <- tapestri.cnv.loh.genotypes(lfile)
  }
  genotypes <- f$genotypes
  kept_cells <- f$cells
  kept_variants <- f$variants
  devnull <- base::gc()

  rownames(genotypes) <- rownames(bc.df)[kept_cells]

  cellcount <- nrow(genotypes)
  varcount <- ncol(genotypes)

  ## Filtering - reduce the genotypic matrix
  ## apply genotype call rate to be present in at least 95% of all cells
  ## exclude if >99% homozygous reference (WT)
  if (z.filter == TRUE) {
    gt.indexes <- tapestri.cnv.loh.filter(genotypes, z.variants, z.missing, z.het.min)
    genotypes.l2 <- genotypes[, gt.indexes]
    kept_variants <- kept_variants[gt.indexes]
  } else {
    genotypes.l2 <- genotypes
  }

  devnull <- base::gc()
  vafmat <- lfile$layers$AD[kept_cells, kept_variants] / lfile$layers$DP[kept_cells, kept_variants]
  devnull <- base::gc()

  rownames(vafmat) <- rownames(genotypes.l2)
  colnames(vafmat) <- colnames(genotypes.l2)

  vafmat.df <- reshape2::melt(vafmat)
  colnames(vafmat.df) <- c("cells", "variants", "af")
  # vafmat.df$variants <- gsub(":", ".", vafmat.df$variants)
  # vafmat.df$variants <- gsub("/", ".", vafmat.df$variants)

  gt.rows <- nrow(genotypes.l2)
  gt.cols <- ncol(genotypes.l2)

  if (gt.rows == 0 | gt.cols == 0) {
    stop("All variants are filtered. Try different filtering settings.")
  }

  ## melt table
  genotypes.fin <- reshape2::melt(as.matrix(genotypes.l2)) # reshape into dataframe
  colnames(genotypes.fin) <- c("cells", "variants", "zygosity")
  genotypes.fin$variants <- gsub(":", ".", genotypes.fin$variants)
  genotypes.fin$variants <- gsub("/", ".", genotypes.fin$variants)

  # merge with amplicons and barcode information from tsv files
  panelannots$varid <- gsub(":", ".", panelannots$varid)
  panelannots$varid <- gsub("/", ".", panelannots$varid)
  genotypes.fin2 <- merge(genotypes.fin, panelannots, by.x = "variants", by.y = "varid")
  genotypes.tab <- merge(genotypes.fin2, vafmat.df, by = c("variants", "cells"))
  # genotypes.tab$cells <- rownames(bc.df)[kept_cells]

  devnull <- base::gc()

  return(genotypes.tab)
}

# Genotype filtering based on cell abundance and allele frquency
#
# @param genotypes genotype data
# @param z.variants zygosity percentage cutoff (Default 90)
# @param z.missing missing zygosity cutoff (Default 100)
# @param z.het.min minimum het percentage (Default 5)
# @return filtered genotypes
# @examples
# \dontrun{
# genotypes.filtered <- tapestri.cnv.loh.filter(genotypes, z.variants, z.missing)
# }
tapestri.cnv.loh.filter <- function(genotypes, z.variants = 90, z.missing = 100, z.het.min = 5) {
  ## Keep columns with variants in more than 5% of the cells
  gt <- list()
  gt$col_homref <- colSums(genotypes == 0)
  gt$col_hetalele <- colSums(genotypes == 1)
  gt$col_homalele <- colSums(genotypes == 2)
  gt$col_nocall <- colSums(genotypes == 3)
  gt.df <- do.call(rbind.data.frame, gt)
  alleles <- colSums(gt.df[1:3, ])
  alleles[alleles == 0] <- 1
  colnames(gt.df) <- colnames(genotypes)
  gt.df.f <- as.matrix(gt.df * 100 / alleles)
  gt.indexes <- gt.df.f[1, ] < z.variants & gt.df.f[2, ] < z.variants & gt.df.f[2, ] > z.het.min & gt.df.f[3, ] < z.variants & gt.df.f[4, ] < z.missing

  return(gt.indexes)
}


# DEPRECATED
# Heatmap for filtered genotype and cells
# 
# @param loh_data loh data
# @param max_cluster max_cluster
# @param outfile output file
# @param freq Plot genotypes with allele frequency
# @param merge_missing Merge missing to WT
# @param median.filter apply median filter while drawing heatmap
# @param median.bin smoothing bin size
# @param filter.common remove variants having similar genotypes between control group and test group
# @param delta.common minimum difference between control/test for variant to be retained
# @return clust
# @export
# @examples
# \dontrun{
# h.clust <- tapestri.cnv.loh.heatmap(loh_data, max_cluster, outfile, TRUE, FALSE)
# }
tapestri.cnv.loh.heatmap <- function(loh_data, max_cluster, outfile, freq = FALSE, merge_missing = FALSE, median.filter = TRUE, median.bin = 11, filter.common = TRUE, delta.common = 0.01) {
  col_name <- "zygosity"
  gt.df <- as.data.frame(loh_data[, c(1, 2, 3)])

  if (freq == TRUE) {
    col_name <- "af"
    gt.df <- as.data.frame(loh_data[, c(1, 2, 13)])
  }

  colnames(gt.df) <- c("variants", "cells", col_name)

  cellcount <- length(unique(gt.df$cells))

  missing_gt <- 3
  legend_params <- list(title = "Genotype", at = c(0, 1, 2, 3), border = "black", labels = c("WT", "HET", "HOM", "Missing"), title_gp = grid::gpar(fontsize = 30, fontface = "bold"), labels_gp = grid::gpar(fontsize = 25), legend_height = 10, legend_width = 20, grid_height = grid::unit(10, "mm"), grid_width = grid::unit(20, "mm"))
  col_fun <- c("grey", "red", "brown", "white")
  if (merge_missing == TRUE) {
    missing_gt <- 0
    legend_params <- list(title = "Genotype", at = c(0, 1, 2), border = "black", labels = c("WT", "HET", "HOM"), title_gp = grid::gpar(fontsize = 30, fontface = "bold"), labels_gp = grid::gpar(fontsize = 25), legend_height = 10, legend_width = 20, grid_height = grid::unit(10, "mm"), grid_width = grid::unit(20, "mm"))
    col_fun <- c("grey", "red", "brown")
  }

  if (freq == FALSE) {
    gt.df$zygosity[gt.df$zygosity == 3] <- missing_gt
  } else {
    gt.df$af[gt.df$af == 1] <- 2
    gt.df$af[gt.df$af < 0.4] <- 0
    gt.df$af[gt.df$af >= 0.4 & gt.df$af < 0.6] <- 1
    gt.df$af[gt.df$af >= 0.6 & gt.df$af < 1] <- 2
    gt.df$af[is.na(gt.df$af)] <- missing_gt
  }
  # gt.df$cells <- c(1:cellcount)

  gt.mat <- as.data.frame(reshape2::acast(gt.df, cells ~ variants, value.var = col_name))
  sorted_col_names <- tapestri.cnv.loh.sort_by_loc(colnames(gt.mat))
  gt.mat <- gt.mat[, sorted_col_names]

  if (filter.common == TRUE) {
    gt.mat <- tapestri.cnv.loh.selectDistinct(gt.mat, max_cluster, delta.common, dirname(outfile))
  }

  sorted_col_names <- colnames(gt.mat)

  distfun <- function(x) stats::dist(x, method = "binary")
  gt.dist <- stats::dist(as.matrix(gt.mat) == 1, "binary")
  gt.clust <- stats::hclust(gt.dist, "average")
  gt.dend <- stats::order.dendrogram(stats::as.dendrogram(gt.clust))

  if (median.filter == TRUE) {
    gt.med <- gt.mat[gt.dend, ]
    gt.med <- apply(gt.med, 2, function(x) stats::runmed(x, median.bin))
    rownames(gt.med) <- rownames(gt.mat)[gt.dend]
    gt.med <- gt.med[rownames(gt.mat), ]
  } else {
    gt.med <- gt.mat
  }

  cell_names <- rownames(gt.med)

  # rownames(gt.med) <- NULL

  var_col_names <- tapestri.cnv.loh.get_column_annots(loh_data, sorted_col_names)
  ## Heatmap visualisation after filtering
  h <- ComplexHeatmap::Heatmap(
    as.matrix(gt.med),
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    row_order = gt.dend,
    column_order = sorted_col_names,
    column_names_centered = TRUE,
    row_km = max_cluster, row_km_repeats = 100,
    row_title = "cluster_%s",
    row_title_gp = grid::gpar(fontsize = 30),
    column_title_gp = grid::gpar(fontsize = 25),
    column_names_gp = grid::gpar(fontsize = 12),
    heatmap_legend_param = legend_params,
    gap = grid::unit(5, "mm"),
    column_gap = grid::unit(5, "mm"),
    use_raster = TRUE,
    raster_quality = 2,
    col = col_fun,
    column_split = factor(gsub("\\..*", "", sorted_col_names), levels = unique(gsub("\\..*", "", sorted_col_names))),
    column_labels = var_col_names,
    cluster_column_slices = FALSE,
    border = TRUE,
    show_row_dend = FALSE,
    show_row_names = FALSE
  )
  grDevices::pdf(outfile, 45, 18)
  h.draw <- ComplexHeatmap::draw(h)
  grDevices::dev.off()
  clust <- snv_clustering(gt.mat, max_cluster)
  return(clust)
}


# Raw genotypes
#
# @param lfile loom data
# @return raw genotypes
# @examples
# \dontrun{
# genotypes <- tapestri.cnv.loh.genotypes(lfile)
# }
tapestri.cnv.loh.genotypes <- function(lfile) {
  allsites_block <- hdf5r::h5attributes(lfile)$allsites_block

  if (length(allsites_block) == 0) {
    genotypes <- lfile$matrix[, ]
    colnames(genotypes) <- as.character(lfile$row.attrs$id[])
  } else {
    genotypes <- lfile$matrix[, 1:allsites_block]
    colnames(genotypes) <- (as.character(lfile$row.attrs$id[]))[1:allsites_block]
  }

  var_ids <- colnames(genotypes)
  var_ids <- gsub(":", ".", var_ids)
  var_ids <- gsub("/", ".", var_ids)
  colnames(genotypes) <- var_ids

  cells <- c(1:dim(genotypes)[1])
  variants <- c(1:dim(genotypes)[2])

  return(list("genotypes" = genotypes, "cells" = cells, "variants" = variants))
}


# Variant annotations
#
# @param lfile loom data
# @return annotations
# @examples
# \dontrun{
# annot <- tapestri.cnv.loh.annotations(lfile)
# }
tapestri.cnv.loh.annotations <- function(lfile) {
  # vectors with information from rows
  panelannots <- data.frame(
    chrom = lfile$row.attrs$CHROM[],
    pos = as.numeric(lfile$row.attrs$POS[]),
    ref = lfile$row.attrs$REF[],
    alt = lfile$row.attrs$ALT[],
    qual = lfile$row.attrs$QUAL[],
    amplicon = as.character(lfile$row.attrs$amplicon[]),
    varid = as.character(lfile$row.attrs$id[]),
    reflen = nchar(lfile$row.attrs$REF[]),
    altlen = nchar(lfile$row.attrs$ALT[]),
    varType = NA
  )

  panelannots$varid <- gsub(":", ".", panelannots$varid)
  panelannots$varid <- gsub("/", ".", panelannots$varid)

  # compute variant type
  panelannots$varType[which(panelannots$reflen == 1 & panelannots$altlen == 1)] <- "SNV"
  panelannots$varType[which(panelannots$reflen == 1 & panelannots$altlen == 0)] <- "MNV"
  panelannots$varType[which(panelannots$reflen > 1 & panelannots$altlen == 1)] <- "DEL"
  panelannots$varType[which(panelannots$reflen == 1 & panelannots$altlen > 1)] <- "INS"
  return(panelannots)
}


# Filter raw genotypes based on quality
#
# @param loom Loom file
# @param gqc Genotype quality cutoff (default 30)
# @param dpc Read depth cutoff (default 10)
# @param afc Allele frequency cutoff (default 20)
# @param mv Remove variants with < mv of known values (default 50)
# @param mc Remove variants with < mc of known values (default 50)
# @param mm Remove variants mutated in < mm of cells (default 1)
# @param gt.mask mask low quality GT as missing (if GQ/DP/AF lower than cutoff, default FALSE)
# @return Filtered genotypes
# @examples
# \dontrun{
# f <- tapestri.cnv.loh.genotypes.filter(loom)
# f <- tapestri.cnv.loh.genotypes.filter(loom, 30, 10, 20, 50, 50, 1)
# }
tapestri.cnv.loh.genotypes.filter <- function(loom, gqc = 30, dpc = 10, afc = 20, mv = 50, mc = 50, mm = 1, gt.mask = FALSE) {
  allsites_block <- hdf5r::h5attributes(loom)$allsites_block

  mask <- (loom$matrix[, 1:allsites_block] < 3)
  devnull <- base::gc()

  gq <- (loom$layers$GQ[, 1:allsites_block] >= gqc)
  devnull <- base::gc()

  gt <- as.matrix(loom$matrix[, 1:allsites_block])
  mutated <- (gt == 1 | gt == 2)
  base::rm(gt)
  devnull <- base::gc()

  dp <- loom$layers$DP[, 1:allsites_block]
  ad <- loom$layers$AD[, 1:allsites_block]
  af <- matrix(100, nrow = nrow(dp), ncol = ncol(dp))
  af[mutated] <- ad[mutated] * 100 / dp[mutated]
  af[is.na(af)] <- 0
  af <- (af >= afc)
  base::rm(ad)
  devnull <- base::gc()
  dp <- (dp >= dpc)
  devnull <- base::gc()

  ngt_filter <- gq & dp & af & mask
  mv.c <- base::colMeans(ngt_filter, na.rm = T) * 100
  kept_variants <- base::which(mv.c >= mv)
  mc.c <- base::rowMeans(ngt_filter[, kept_variants], na.rm = T) * 100
  kept_cells <- base::which(mc.c >= mc)

  ngt_mutated <- mutated & ngt_filter

  ngt_mutated <- ngt_mutated[kept_cells, kept_variants]
  mm.c <- base::colMeans(ngt_mutated, na.rm = T) * 100

  mv.f <- (mv.c < mv)
  mm.f <- (mm.c < mm)
  mv.f[!mv.f] <- mm.f

  kept_variants <- base::which(!mv.f)
  if (!(length(kept_variants) & length(kept_cells))) {
    stop("All cells/variants are filtered. Try different filtering settings.")
  }

  gt <- loom$matrix[, 1:allsites_block]
  if (gt.mask == TRUE) {
    gt[!ngt_filter] <- 3
  }
  gt <- gt[kept_cells, kept_variants]
  var_ids <- (as.character(loom$row.attrs$id[]))[kept_variants]
  var_ids <- gsub(":", ".", var_ids)
  var_ids <- gsub("/", ".", var_ids)
  colnames(gt) <- var_ids

  devnull <- base::gc()
  return(list("genotypes" = gt, "cells" = kept_cells, "variants" = kept_variants))
}

# Sort column names by chromosome and location
#
# @param col_names column names
# @return Sorted column names
# @examples
# \dontrun{
# col_names <- tapestri.cnv.loh.sort_by_loc(colnames(genotypes.df))
# }
tapestri.cnv.loh.sort_by_loc <- function(col_names) {
  out <- strsplit(as.character(col_names), '\\.')
  out <- do.call(rbind, out)
  out[out == "chrX"] <- "chr23"
  out[out == "chrY"] <- "chr24"
  out[, 1] <- gsub("chr", "", out[, 1])
  out <- as.data.frame(out)
  out$V1 <- as.numeric(as.character(out$V1))
  out$V2 <- as.numeric(as.character(out$V2))
  out <- out[order(out$V1, out$V2), ]
  out$V1 <- as.character(out$V1)
  out$V2 <- as.character(out$V2)
  out$V1[out$V1 == 23] <- "X"
  out$V1[out$V1 == 24] <- "Y"
  col_names <- paste(out$V1, out$V2, out$V3, out$V4, sep = ".")
  col_names <- paste("chr", col_names, sep = "")
  return(col_names)
}

# Select Distinctive Mutations
#
# @param gt.mat filtered data
# @param max_cluster number of cluster to be evaluated
# @param delta.common difference to allow variants
# @param outdir directory to store temporary heatmap output
# @return gt.mat
# @examples
# \dontrun{
# gt.mat <- tapestri.cnv.loh.selectDistinct(gt.mat)
# }
tapestri.cnv.loh.selectDistinct <- function(gt.mat, max_cluster = 2, delta.common = 0.01, outdir = NULL) {
  h <- ComplexHeatmap::Heatmap(as.matrix(gt.mat), row_km = max_cluster, row_km_repeats = 100)

  # need a better logic to disable plotting
  temp_file <- tempfile(pattern = "Clust_", tmpdir = outdir, fileext = ".pdf")
  grDevices::pdf(temp_file)
  h.draw <- ComplexHeatmap::draw(h)
  grDevices::dev.off()
  if (file.exists(temp_file)) {
    file.remove(temp_file)
  }

  h.clust <- ComplexHeatmap::row_order(h.draw)

  het_sums <- rep.int(0, length(h.clust))
  for (i in c(1:length(h.clust))) {
    het_sums[i] <- sum(gt.mat[h.clust[[i]], ] == 1)
  }
  norm.clust <- order(het_sums, decreasing = TRUE)[1]
  test.clusts <- order(het_sums, decreasing = TRUE)[-1]

  control.df <- gt.mat[unlist(h.clust[norm.clust]), ]
  test.df <- gt.mat[unlist(h.clust[-norm.clust]), ]

  control.mean <- apply(control.df, 2, mean)
  test.mean <- apply(test.df, 2, mean)

  delta <- abs(control.mean - test.mean)

  gt.mat <- gt.mat[, delta > delta.common]

  return(gt.mat)
}


# Annotate variants with genes
#
# @param loh loh data
# @param sorted_col_names sorted variant column names
# @return annotated variant names
# @examples
# \dontrun{
# variant_annots <- tapestri.cnv.loh.get_column_annots(loh, sorted_col_names)
# }
tapestri.cnv.loh.get_column_annots <- function(loh, sorted_col_names) {
  variant_regions <- unique(dplyr::select(loh, "chrom", "pos", "variants"))
  colnames(variant_regions)[3] <- "name"
  variant_regions$start <- variant_regions$pos - 1
  colnames(variant_regions)[2] <- "end"
  variant_regions <- variant_regions[, c("chrom", "start", "end", "name")]
  variant_regions$chrom <- paste("chr", variant_regions$chrom, sep = "")
  variant_regions <- regioneR::toGRanges(variant_regions)

  # annots <- c('hg19_basicgenes')
  # annotations <- annotatr::build_annotations(genome = 'hg19', annotations = 'hg19_basicgenes')
  # dm_annotated <- annotatr::annotate_regions(regions = variant_regions, annotations = annotations, ignore.strand = TRUE, quiet = TRUE)
  # dm_annotated <- data.frame(dm_annotated)
  # dm_annotated <- dplyr::select(dm_annotated, "annot.symbol", "name")
  # dm_annotated <- unique(dm_annotated)
  # dm_annotated <- dm_annotated[!is.na(dm_annotated$annot.symbol),]
  # ordered.variants <- as.data.frame(sorted_col_names)
  # colnames(ordered.variants) <- "id"
  # ordered.variants$id <- as.character(ordered.variants$id)
  # ordered.variants <- base::merge(ordered.variants, dm_annotated, by.x="id", by.y="variants", all.x=TRUE)
  # ordered.variants <- ordered.variants[match(sorted_col_names, ordered.variants$id),]
  # variant_annots <- ifelse(is.na(ordered.variants$annot.symbol), ordered.variants$id, ordered.variants$annot.symbol)
  ordered.variants <- annotate_regions_with_genes(variant_regions, sorted_col_names)
  variant_annots <- ifelse(is.na(ordered.variants$annot.symbol), ordered.variants$id, ordered.variants$annot.symbol)
  return(variant_annots)
}

#' Get cell cluster labels
#' @param genotypes filtered gt matrix
#' @param max_cluster maximum number of clusters
#' @return cell cluster labels as dataframe
#' @export
#' @examples
#' \dontrun{
#' snv_clust <- snv_clustering(genotypes, max_cluster)
#' }
snv_clustering <- function(genotypes, max_cluster) {
  set.seed(42)

  sorted_col_names <- tapestri.cnv.loh.sort_by_loc(colnames(genotypes))
  genotypes <- genotypes[, sorted_col_names]

  clusters <- stats::kmeans(genotypes, max_cluster, iter.max = 100)
  clus <- as.data.frame(as.factor(clusters$cluster))
  colnames(clus) <- "cluster"
  return(clus)
}


#' Extract genotypes from loom
#'
#' @param lfile loom file
#' @param barcodes normalized barcode distribution data
#' @param gt.filter enable genotype quality filter (Enabled by default)
#' @param gt.gqc Genotype quality cutoff (default 30)
#' @param gt.dpc Read depth cutoff (default 10)
#' @param gt.afc Allele frequency cutoff (default 20)
#' @param gt.mv Remove variants genotyped in < X percent of cells (default 50)
#' @param gt.mc Remove cells with mutations in < X percent of variants (default 50)
#' @param gt.mm Remove variants mutated in < X percent of cells (default 1)
#' @param gt.mask mask low quality GT as missing (if GQ/DP/AF lower than cutoff, default FALSE)
#' @return loh data
#' @export
#' @examples
#' \dontrun{
#' genotypes <- extract_genotypes(lfile, barcodes, gt.filter = FALSE)
#' genotypes <- extract_genotypes(lfile, barcodes,
#'   gt.filter = TRUE, gt.gqc = 30, gt.dpc = 10, gt.afc = 20,
#'   gt.mv = 50, gt.mc = 50, gt.mm = 1, gt.mask = FALSE
#' )
#' }
extract_genotypes <- function(lfile, barcodes, gt.filter = TRUE, gt.gqc = 30, gt.dpc = 10, gt.afc = 20, gt.mv = 50, gt.mc = 50, gt.mm = 1, gt.mask = FALSE) {

  # matrix with genotypic information
  if (gt.filter == TRUE) {
    f <- tapestri.cnv.loh.genotypes.filter(lfile, gqc = gt.gqc, dpc = gt.dpc, afc = gt.afc, mv = gt.mv, mc = gt.mc, mm = gt.mm, gt.mask = gt.mask)
  } else {
    f <- tapestri.cnv.loh.genotypes(lfile)
  }
  genotypes <- f$genotypes
  kept_cells <- f$cells
  kept_variants <- f$variants
  devnull <- base::gc()

  rownames(genotypes) <- rownames(barcodes)[kept_cells]
  return(genotypes)
}

#' Filter genotypes based on zygpsity fraction
#'
#' @param genotypes genotype data
#' @param z.variants zygosity percentage cutoff (Default 90)
#' @param z.missing missing zygosity cutoff (Default 100)
#' @param z.het.min minimum het percentage (Default 5)
#' @return filtered genotypes
#' @export
#' @examples
#' \dontrun{
#' genotypes.filtered <- filter_by_zygosity(genotypes, 90, 100, 5)
#' }
filter_by_zygosity <- function(genotypes, z.variants = 90, z.missing = 100, z.het.min = 5) {
  ## Keep columns with variants in more than 5% of the cells
  gt <- list()
  gt$col_homref <- colSums(genotypes == 0)
  gt$col_hetalele <- colSums(genotypes == 1)
  gt$col_homalele <- colSums(genotypes == 2)
  gt$col_nocall <- colSums(genotypes == 3)
  gt.df <- do.call(rbind.data.frame, gt)
  alleles <- colSums(gt.df[1:3, ])
  alleles[alleles == 0] <- 1
  colnames(gt.df) <- colnames(genotypes)
  gt.df.f <- as.matrix(gt.df * 100 / alleles)
  gt.indexes <- gt.df.f[1, ] < z.variants & gt.df.f[2, ] < z.variants & gt.df.f[2, ] > z.het.min & gt.df.f[3, ] < z.variants & gt.df.f[4, ] < z.missing

  genotypes <- genotypes[, gt.indexes]
  return(genotypes)
}

#' Select Distinctive Mutations
#'
#' @param genotypes filtered data
#' @param max_cluster number of cluster to be evaluated
#' @param delta.common difference to allow variants
#' @return filtered genotypes
#' @export
#' @examples
#' \dontrun{
#' genotypes.f <- filter_common_variants(genotypes)
#' }
filter_common_variants <- function(genotypes, max_cluster = 2, delta.common = 0.01) {
  h <- ComplexHeatmap::Heatmap(as.matrix(genotypes), row_km = max_cluster, row_km_repeats = 100)

  # need a better logic to disable plotting
  temp_file <- tempfile(pattern = "Clust_", tmpdir = ".", fileext = ".pdf")
  grDevices::pdf(temp_file)
  h.draw <- ComplexHeatmap::draw(h)
  grDevices::dev.off()
  if (file.exists(temp_file)) {
    file.remove(temp_file)
  }

  h.clust <- ComplexHeatmap::row_order(h.draw)

  het_sums <- rep.int(0, length(h.clust))
  for (i in c(1:length(h.clust))) {
    het_sums[i] <- sum(genotypes[h.clust[[i]], ] == 1)
  }
  norm.clust <- order(het_sums, decreasing = TRUE)[1]
  test.clusts <- order(het_sums, decreasing = TRUE)[-1]

  control.df <- genotypes[unlist(h.clust[norm.clust]), ]
  test.df <- genotypes[unlist(h.clust[-norm.clust]), ]

  control.mean <- apply(control.df, 2, mean)
  test.mean <- apply(test.df, 2, mean)

  delta <- abs(control.mean - test.mean)

  genotypes <- genotypes[, delta > delta.common]

  return(genotypes)
}

#' Annotate genotypes
#'
#' @param loom loom file object
#' @param barcodes barcodes data
#' @param genotypes genotype data
#' @return annotated genotypes
#' @export
#' @examples
#' \dontrun{
#' genotypes.annotated <- annotate_genotypes(loom, barcodes, genotypes)
#' }
annotate_genotypes <- function(loom, barcodes, genotypes) {
  kept_cells <- match(rownames(genotypes), rownames(barcodes))
  var_ids <- (as.character(loom$row.attrs$id[]))
  var_ids <- gsub(":", ".", var_ids)
  var_ids <- gsub("/", ".", var_ids)
  kept_variants <- match(colnames(genotypes), var_ids)

  vafmat <- loom$layers$AD[kept_cells, kept_variants] / loom$layers$DP[kept_cells, kept_variants]
  devnull <- base::gc()

  rownames(vafmat) <- rownames(genotypes)
  colnames(vafmat) <- colnames(genotypes)

  vafmat.df <- reshape2::melt(vafmat)
  colnames(vafmat.df) <- c("cells", "variants", "af")

  genotypes.df <- reshape2::melt(as.matrix(genotypes))

  colnames(genotypes.df) <- c("cells", "variants", "zygosity")

  # merge with amplicons and barcode information from tsv files
  panelannots <- tapestri.cnv.loh.annotations(loom)
  genotypes.df <- merge(genotypes.df, panelannots, by.x = "variants", by.y = "varid")
  genotypes.df <- merge(genotypes.df, vafmat.df, by = c("variants", "cells"))

  devnull <- base::gc()
  return(genotypes.df)
}

#' Heatmap for filtered genotype and cells
#'
#' @param genotypes genotypes
#' @param max_cluster max_cluster
#' @param outfile output file
#' @param median.filter apply median filter while drawing heatmap
#' @param median.bin smoothing bin size
#' @param annotations annotated genotypes dataframe
#' @return clust
#' @export
#' @examples
#' \dontrun{
#' clust <- snv_heatmap(genotypes, max_cluster, outfile, TRUE, FALSE)
#' }
snv_heatmap <- function(genotypes, max_cluster=NULL, outfile=NULL, median.filter = TRUE, median.bin = 11, annotations = NULL) {
  if (is.null(max_cluster)) {
    stop("max_cluster not provided")
  }
  legend_params <- list(title = "Genotype", at = c(0, 1, 2, 3), border = "black", labels = c("WT", "HET", "HOM", "Missing"), title_gp = grid::gpar(fontsize = 8, fontface = "bold"), labels_gp = grid::gpar(fontsize = 6), legend_height = 4, legend_width = 6, grid_height = grid::unit(4, "mm"), grid_width = grid::unit(8, "mm"))
  col_fun <- c("grey", "red", "brown", "white")

  sorted_col_names <- tapestri.cnv.loh.sort_by_loc(colnames(genotypes))
  genotypes <- genotypes[, sorted_col_names]

  sorted_col_names <- colnames(genotypes)

  distfun <- function(x) stats::dist(x, method = "binary")
  gt.dist <- stats::dist(as.matrix(genotypes) == 1, "binary")
  gt.clust <- stats::hclust(gt.dist, "average")
  gt.dend <- stats::order.dendrogram(stats::as.dendrogram(gt.clust))

  if (median.filter == TRUE) {
    gt.med <- genotypes[gt.dend, ]
    gt.med <- apply(gt.med, 2, function(x) stats::runmed(x, median.bin))
    rownames(gt.med) <- rownames(genotypes)[gt.dend]
    gt.med <- gt.med[rownames(genotypes), ]
  } else {
    gt.med <- genotypes
  }

  cell_names <- rownames(gt.med)

  var_col_names <- sorted_col_names
  if (!is.null(annotations)) {
    var_col_names <- tapestri.cnv.loh.get_column_annots(annotations, sorted_col_names)
  }
  ## Heatmap visualisation after filtering
  h <- ComplexHeatmap::Heatmap(
    as.matrix(gt.med),
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    row_order = gt.dend,
    column_order = sorted_col_names,
    column_names_centered = TRUE,
    row_km = max_cluster, row_km_repeats = 100,
    row_title = "cluster_%s",
    row_title_gp = grid::gpar(fontsize = 10),
    column_title_gp = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 4),
    heatmap_legend_param = legend_params,
    gap = grid::unit(2, "mm"),
    column_gap = grid::unit(2, "mm"),
    use_raster = TRUE,
    raster_quality = 2,
    col = col_fun,
    column_split = factor(gsub("\\..*", "", sorted_col_names), levels = unique(gsub("\\..*", "", sorted_col_names))),
    column_labels = var_col_names,
    cluster_column_slices = FALSE,
    border = TRUE,
    show_row_dend = FALSE,
    show_row_names = FALSE,
    heatmap_width = grid::unit(30, "cm"),
    heatmap_height = grid::unit(15, "cm")
  )
  if (!is.null(outfile)) {
    grDevices::pdf(outfile, 15, 7)
  }
  h.draw <- ComplexHeatmap::draw(h)

  if (!is.null(outfile)) {
    grDevices::dev.off()
  }
  return(h)
}

#' Genotypes from allele frequency
#'
#' @param loh Annotated genotypes
#' @return annotated genotypes
#' @export
#' @examples
#' \dontrun{
#' genotypes <- genotypes_from_allele_frequency(loh)
#' }
genotypes_from_allele_frequency <- function(loh) {
  missing_gt <- 3
  gt.df <- loh[, c("variants", "cells", "af")]
  gt.df$af[gt.df$af == 1] <- 2
  gt.df$af[gt.df$af < 0.4] <- 0
  gt.df$af[gt.df$af >= 0.4 & gt.df$af < 0.6] <- 1
  gt.df$af[gt.df$af >= 0.6 & gt.df$af < 1] <- 2
  gt.df$af[is.na(gt.df$af)] <- missing_gt
  genotypes <- as.matrix(reshape2::acast(gt.df, cells ~ variants, value.var = "af"))
  return(genotypes)
}
