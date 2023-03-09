#' Connect to LOOM file
#'
#' @param loom_file LOOM file
#' @return loom connection to input LOOM file
#' @export
#' @examples
#' \dontrun{
#' loom <- connect_to_loom(loom_file)
#' }
connect_to_loom <- function(loom_file) {
  loom <- loomR::connect(filename = loom_file, mode = "r+")
  return(loom)
}
