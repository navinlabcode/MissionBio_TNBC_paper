#' Load bed data
#'
#' @param bedfile bed file
#' @return bed records
#' @export
#' @examples
#' \dontrun{
#' ampl <- parse_bed("Myeloid_v2_amplicon.bed")
#' targets <- parse_bed("Myeloid_v2_submitted.bed")
#' }
parse_bed <- function(bedfile) {
  bed <- utils::read.table(paste0(bedfile), sep = "\t", header = F, skip = 1)
  colnames(bed) <- c("chr", "start", "end", "name")
  return(bed)
}
