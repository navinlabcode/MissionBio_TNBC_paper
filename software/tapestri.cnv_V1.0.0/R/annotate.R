# Annotate region list
# @param regions region dataframe from parseBed
# @param sorted_col_names column names soorted by chromosome position
# @return annotated column names
# @examples
# \dontrun{
# ordered.variants <- annotate_regions_with_genes(regions, sorted_col_names)
# }

annotate_regions_with_genes <- function(regions, sorted_col_names) {
  annots <- c("hg19_basicgenes")
  annotations <- annotatr::build_annotations(genome = "hg19", annotations = "hg19_basicgenes")
  dm_annotated <- annotatr::annotate_regions(regions = regions, annotations = annotations, ignore.strand = TRUE, quiet = TRUE)
  dm_annotated <- data.frame(dm_annotated)
  dm_annotated <- dplyr::select(dm_annotated, "annot.symbol", "name")
  dm_annotated <- unique(dm_annotated)
  dm_annotated <- dm_annotated[!is.na(dm_annotated$annot.symbol), ]
  ordered.variants <- as.data.frame(sorted_col_names)
  colnames(ordered.variants) <- "id"
  ordered.variants$id <- as.character(ordered.variants$id)
  ordered.variants <- base::merge(ordered.variants, dm_annotated, by.x = "id", by.y = "name", all.x = TRUE)
  ordered.variants <- ordered.variants[match(sorted_col_names, ordered.variants$id), ]
  # col_annots <- ifelse(is.na(ordered.variants$annot.symbol), ordered.variants$id, ordered.variants$annot.symbol)
  return(ordered.variants)
}
