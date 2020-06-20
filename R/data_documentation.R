#' Example of object obtained from using the master_matrix function
#'
#' A list of 4 elements (m_matrix, polygon, raster_base, and PCA_results).
#' See function \code{\link{master_matrix}}.
#'
#' @format A list of 4 elements:
#' \describe{
#'   \item{m_matrix}{data frame wiht 6276 rows and 10 columns}
#'   \item{polygon}{object of class SpatialPolygonsDataFrame}
#'   \item{raster_base}{object of class RasterLayer}
#'   \item{PCA_results}{list of length 5}
#' }
#'
#' @examples
#' data("m_matrix", package = "biosurvey")
#'
#' summary(m_matrix)
"m_matrix"


#' Example of spatial polygon for a region of interest
#'
#' An object of class SpatialPolygonsDataFrame.
#'
#' @format SpatialPolygonsDataFrame:
#' \describe{
#'   \item{data}{data frame wiht 1 row and 11 columns}
#'   \item{polygons}{SpatialPolygons}
#'   \item{proj4string}{object of class CRS}
#' }
#'
#' @examples
#' data("mx", package = "biosurvey")
#'
#' mx
"mx"


#' Example of a data frame of species' found in distint positions
#'
#' @description A data frame with two columns: "ID" and "Species" and 590 rows.
#' It is the output of the function \code{\link{spdf_2data}}.
#'
#' @format data frame:
#' \describe{
#'   \item{ID}{identifier of position}
#'   \item{Species}{different species in the table}
#' }
#'
#' @examples
#' data("sp_data", package = "biosurvey")
#'
#' summary(sp_data)
"sp_data"


#' Example of species ranges as SpatialPolygonsDataFrame
#'
#' An object of class SpatialPolygonsDataFrame.
#'
#' @format data frame:
#' \describe{
#'   \item{data}{data frame wiht 25 rows and 1 column}
#'   \item{polygons}{SpatialPolygons}
#'   \item{proj4string}{object of class CRS}
#' }
#'
#' @examples
#' data("species_data", package = "biosurvey")
#'
#' species_data
"species_data"


#' Example of stack of layers of suitable and unsuitable conditions for species
#'
#' An object of class RasterStack containing information about suitable and
#' unsuitable conditions for five species.
#'
#' @format A RasterStack with 109 rows, 182 columns, 19838 cells, and 5 layers:
#' \describe{
#'   \item{RasterLayer}{suitable (1) and unsuitable (0) conditions}
#' }
#'
#' @examples
#' sp_layers <- raster::stack(system.file("extdata/sp_layers.tif",
#'                            package = "biosurvey"))
#'
#' sp_layers
#' @name sp_layers
NULL


#' Example of variables to be used for preparing a master matrix
#'
#' A dataset containing raster variables for an area that is relevant for used
#' in examples included in the package \code{\link{biosurvey}}.
#'
#' @format A RasterStack with 109 rows, 190 columns, 20710 cells, and 6 layers:
#' \describe{
#'   \item{Temperature}{temperature, in Celsius degrees times 10.}
#'   \item{Precipitation}{precipitation, in milimeters.}
#' }
#'
#' @source \url{https://www.worldclim.org/data/index.html}
#'
#' @examples
#' variables <- raster::stack(system.file("extdata/variables.tif",
#'                                        package = "biosurvey"))
#'
#' raster::plot(variables[[1]])
#' @name variables
NULL


#' Occurrence records for the species Parides gundlachianus
#'
#' @description A dataset containing geographic coordinates of a Cuban butterflie.
#'
#' @format A data frame with 19 rows and 3 columns.
#' \describe{
#'   \item{name}{character, species scientific name.}
#'   \item{longitude}{numeric, longitude values.}
#'   \item{latitude}{numeric, latitude values.}
#' }
#' @source \url{https://www.gbif.org/}
#'
#' @examples
#' data("sp_occurrences", package = "biosurvey")
#' head(sp_occurrences)
"sp_occurrences"


#' A list of vectors of distances
#'
#' @description A list of six vectors of point distances.
#'
#' @format A list.
#' \describe{
#'   \item{vector}{numeric, values of distances (six elements)}
#' }
#' @source \url{https://www.gbif.org/}
#'
#' @examples
#' data("dist_list", package = "biosurvey")
#' "dist_list"
#'
