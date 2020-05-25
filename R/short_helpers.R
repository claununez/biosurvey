#' Test wheather a number is pair
#'
#' @param x (numeric) value to be tested.
#' @return Logical value
#' @export
#' @examples
#' is_pair(4)
#' is_pair(5)

is_pair <- function(x) {
  x / 2 == as.integer(x / 2)
}





#' Helper function to find raster extension
#'
#' @param format (character) any of the format types allowed for raster objects.
#' See \code{\link[raster]{writeFormats}} (e.g., "GTiff").
#'
#' @return Raster extension according to format type.
#'
#' @export
#'
#' @examples
#' match_rformat("GTiff")

match_rformat <- function(format) {
  if (missing(format)) {stop("Argument 'format needs to be defined.")}
  if (format == "raster") {format1 <- ".grd"}
  if (format == "GTiff") {format1 <- ".tif"}
  if (format == "EHdr") {format1 <- ".bil"}
  if (format == "ascii") {format1 <- ".asc"}
  if (format == "SAGA") {format1 <- ".sdat"}
  if (format == "IDRISI") {format1 <- ".rst"}
  if (format == "CDF") {format1 <- ".nc"}
  if (format == "ENVI") {format1 <- ".envi"}
  if (format == "HFA") {format1 <- ".img"}
  return(format1)
}
