#' Read barcode distribution data
#'
#' @param files List of tab-delimited barcode distribution files output by Tapestri Pipeline
#' @param headerfile vcf header file with barcode mapping
#' @return barcodes - data frame containing merged barcode data from all barcode files
#' @export
#' @examples
#' \dontrun{
#' barcodes <- read_barcodes(
#'   c(
#'     "/path/to/tube-0.barcode.cell.distribution.tsv",
#'     "/path/to/tube-1.barcode.cell.distribution.tsv"
#'   )
#' )
#' }
read_barcodes <- function(files, headerfile) {
  barcodes <- matrix(, nrow=0, ncol=0)

  files <- sort(files)

  for (f in files) {
    bc <- parse_barcode_file(f)
    if (ncol(barcodes) == 0) {
      barcodes <- matrix(, nrow=0, ncol=ncol(bc))
    }
    barcodes <- rbind(as.matrix(barcodes), as.matrix(bc))
  }
  barcodes <- as.data.frame(barcodes)
  bclist <- read.table(headerfile)
  bclist <- t(bclist)
  if(length(bclist) != nrow(barcodes)) {
    stop("Barcodes do not match with vcf header.")
  }
  if (!any(stringr::str_detect(rownames(barcodes), "-\\d$"))) {
    bclist <- gsub("-.*", "", bclist) # remove trailing -1
    rownames(bclist) <- bclist[,1]
    bclist <- as.data.frame(bclist)
    bclist <- rownames(bclist)
  }
  barcodes <- barcodes[match(bclist, rownames(barcodes)), ]
  if (sum(is.na(barcodes)) > 0) {
    stop("Barcodes distribution file has missing data. Please verify distribution file has same barcodes as vcf_header.")
  }
  return(barcodes)
}


#' Normalize barcode data
#'
#' @param barcodes data frame containing merged barcodes from all tubes
#' @param advanced advance normalization (dividing read counts by cell median before applying normalization by total sum)
#' @return normalized_barcodes - normalized barcode data
#' @export
#' @examples
#' \dontrun{
#' normalized_barcodes <- normalize_barcodes(barcodes)
#' }
normalize_barcodes <- function(barcodes, advanced = TRUE) {
  normalized_barcodes <- barcodes

  if (advanced == TRUE) {
    df <- as.matrix(normalized_barcodes)

    y <- colSums(df)
    z <- y / stats::median(y)
    df <- df[, names(z[z > 0.2])]

    ampnorm <- c()

    for (col in colnames(df)) {
      ampnorm <- cbind(ampnorm, as.numeric(round(df[, col] / z[col], 3)))
    }
    colnames(ampnorm) <- colnames(df)
    rownames(ampnorm) <- rownames(df)

    normalized_barcodes <- ampnorm
  }
  normalized_barcodes <- t(apply(
    normalized_barcodes, 1,
    norm <- function(x) {
      return(x / sum(x))
    }
  ))

  return(normalized_barcodes)
}


# Parse barcode distribution file
#
# @param barcode_file Tab-delimited barcode distribution file
# @return barcodes - data frame containing barcode data from input file
# @examples
# \dontrun{
# barcodes <- parse_barcode_file("/path/to/tube-0.barcode.cell.distribution.tsv")
# }
parse_barcode_file <- function(barcode_file) {
  barcodes <- utils::read.table(barcode_file, sep = "\t", header = T, row.names = 1)

  # Some of the barcode distribution files have extra empty columns at the end.
  # This caused the column header to shift incorrectly to right.
  # The code below handles such scenarios
  if (all(is.na(barcodes[, ncol(barcodes)])) && colnames(barcodes)[1] == "X") {
    cnames <- colnames(barcodes)
    barcodes <- barcodes[, 1:ncol(barcodes) - 1]
    colnames(barcodes) <- cnames[2:length(cnames)]
  }
  return(barcodes)
}
