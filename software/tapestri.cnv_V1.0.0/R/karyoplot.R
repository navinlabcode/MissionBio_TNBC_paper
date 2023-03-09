if (getRversion() >= "2.15.1") utils::globalVariables(c("name", "annot.symbol", "ploidy.to.genes.ranges"))

# Plot ploidy ove karyplot
#
# @param ploidy ploidy
# @param amplicons amplicon bed file
# @param outfile outfile
# @param clust clustering output from loh heatmap
# @param cluster_id normal or control cluster
# @param min.amplicons filter genes with amplicons < min.amplicons
# @importFrom magrittr %>%
# @examples
# \dontrun{
# ploidy_karyoplot(ploidy, amplicons, 'karyoplot.pdf')
# }
ploidy_karyoplot <- function(ploidy, amplicons, outfile, clust, cluster_id = NULL, min.amplicons = 3) {
  amplicons <- amplicons[amplicons$name %in% colnames(ploidy), ]
  ploidy.melt <- reshape2::melt(ceiling(as.matrix(ploidy)))
  colnames(ploidy.melt) <- c("cells", "amplicon", "cn")

  if (!is.null(cluster_id)) {
    group.cells <- rownames(clust)[clust$cluster == cluster_id]
    group.ploidy <- ploidy[group.cells, ]
    ploidy.melt <- reshape2::melt(as.matrix(group.ploidy))
    colnames(ploidy.melt) <- c("cells", "amplicon", "cn")
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

  ploidy_to_genes <- ploidy_to_genes[!is.na(ploidy_to_genes$annot.symbol), ]
  ploidy_to_genes <- ploidy_to_genes %>%
    dplyr::group_by(cells, annot.symbol) %>%
    dplyr::filter(dplyr::n_distinct(name) >= min.amplicons)
  ploidy_to_genes <- ploidy_to_genes %>%
    dplyr::group_by(annot.symbol) %>%
    dplyr::mutate("cn_median" = stats::median(ploidy_to_genes$cn, na.rm = T))

  ploidy_to_genes <- as.data.frame(ploidy_to_genes)
  ploidy_to_genes <- unique(ploidy_to_genes[, c("seqnames", "annot.start", "annot.end", "annot.symbol", "cn_median")])
  colnames(ploidy_to_genes) <- c("seqnames", "start", "end", "gene", "cn_median")

  ploidy_to_genes.ranges <- regioneR::toGRanges(ploidy_to_genes)
  ploidy_to_genes.ranges_unique <- GenomicRanges::reduce(ploidy_to_genes.ranges, with.revmap = TRUE, ignore.strand = TRUE)
  for (i in c(1:length(ploidy_to_genes.ranges_unique$revmap))) {
    ploidy_to_genes.ranges_unique$gene[i] <- unique(ploidy_to_genes.ranges$gene[unlist(ploidy_to_genes.ranges_unique$revmap[i])])[1]
    ploidy_to_genes.ranges_unique$cn_median[i] <- unique(ploidy_to_genes.ranges$cn_median[unlist(ploidy_to_genes.ranges_unique$revmap[i])])[1]
  }

  # GenomeInfoDb::seqlevels(ploidy_to_genes.ranges) <- paste0("chr", c(1:22, "X", "Y"))
  ploidy_to_genes.ranges <- GenomeInfoDb::sortSeqlevels(ploidy_to_genes.ranges_unique)
  ploidy_to_genes.ranges <- sort(ploidy_to_genes.ranges)
  ploidy_to_genes.ranges$cn_median <- as.integer(ploidy_to_genes.ranges$cn_median)

  grDevices::pdf(outfile)

  ## plot
  kp <- karyoploteR::plotKaryotype()
  karyoploteR::kpAddBaseNumbers(kp)
  cn_col <- ifelse(ploidy_to_genes.ranges$cn_median > 2, "red", ifelse(ploidy_to_genes.ranges$cn_median < 2, "blue", "black"))
  karyoploteR::kpPlotMarkers(
    kp,
    data = ploidy_to_genes.ranges,
    labels = ploidy_to_genes.ranges$gene,
    label.color = cn_col,
    r1 = 0.7,
    cex = 0.5,
    marker.parts = c(0.2, 0.7, 0.1)
  )
  grDevices::dev.off()
}
