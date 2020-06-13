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





#' Project spatial points from geographic coordinates
#'
#' @param data a matrix or a data frame that contains at least two columns, one with
#' longitude information and the other with latitud information.
#' @param longitude (character) the name of the column that contains the longitude
#' information.
#' @param latitude (character) the name of the column that contains the latitude
#' information.
#' @param which (character) type of projection. There are two options available:
#' "ED", for Azimuthal Equidistant and "EA", for Lambert Azimuthal Equal-Area.
#' Default = "ED".
#'
#' @return SpatialPointsDataFrame projected to an option in \code{which}.
#'
#' @usage
#' wgs84_2aed_laea(data, longitude, latitude, which = "ED")
#'
#' @export
#' @importFrom sp CRS SpatialPointsDataFrame spTransform
#'
#' @examples
#' data("sp_occurrences", package = "biosurvey")
#'
#' sp_occ <- wgs84_2aed_laea(sp_occurrences, longitude = "longitude",
#'                           latitude = "latitude", which = "EA")

wgs84_2aed_laea <- function (data, longitude, latitude, which = "ED") {
  # Initial tests
  if (missing(data)) {
    stop("Argument 'data' is not defined.")
  }
  if (missing(longitude)) {
    stop("Argument 'longitude' is not defined.")
  }
  if (missing(latitude)) {
    stop("Argument 'latitude' is not defined.")
  }

  # Initial projection
  WGS84 <- sp::CRS("+init=epsg:4326")
  dat_s <- sp::SpatialPointsDataFrame(data[, c(longitude, latitude)],
                                      data, proj4string = WGS84)

  # Projecting points
  cent <- apply(data[, c(longitude, latitude)], 2, mean)
  ini <- ifelse(which[1] == "ED", "+proj=aeqd", "+proj=laea")
  prj <- sp::CRS(paste0(ini, " +lat_0=", cent[2], " +lon_0=", cent[1],
                        " +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

  dat_s <- sp::spTransform(dat_s, prj)

  return(dat_s)
}


# Create a bar legend for use in plotting functions

bar_legend <- function (value_range, col, alpha = 1, title = NULL, round = 0) {
  # Initial tests
  if (missing(value_range)) {
    stop("Argument 'value_range' is required to produce the legend.")
  }
  if (missing(col)) {
    stop("Argument 'col' is required to produce the legend.")
  }

  legend_image <- as.raster(matrix(scales::alpha(rev(col), alpha), ncol = 1))
  text(x = 0.6, y = 0.525, labels = title, srt = 90)
  if (is.numeric(value_range)) {
    vals <- round(value_range, round)
  } else {
    vals <- value_range
  }
  text(x = 0.7, y = seq(0.2, 0.85, l = 2), labels = vals, cex = 0.8)
  rasterImage(legend_image, 0.1, 0.2, 0.3, 0.85)
}
