#' Prepare a base object to perform further analyses
#'
#' @description Prepare an S3 object that will serve as the base to perform all
#' further analyses. This object will contain geographic and environmental
#' information that will be used to characterize the region of interest.
#'
#' @param region SpatialPolygonsDataFrame of the region of interest.
#' @param variables RasterStack or RasterBrick of environmental variables.
#' @param do_pca (logical) whether or not to perfor a Principal Component Analysis.
#' Default = FALSE.
#' @param center (logical) whether or not to center variables. Argument to be passed
#' to the function \code{\link[stats]{prcomp}}. Default = TRUE.
#' @param scale (logical) whether or not to scale the variables. Recommended when
#' variables are in different units. Argument to be passed to the function
#' \code{\link[stats]{prcomp}}. Default = FALSE.
#' @param variables_in_matrix (character) name of variables to include in matrix.
#' If NULL (the default) all variables will be inluded.
#'
#' @return
#' A list containing the following elements:
#' - A master matrix for the region of interest, with information
#' about geographic location of raster cells, initial environmental data,
#' and if \code{do_pca} is TRUE, the principal components derived from original data.
#' - A spatial polygon representing the region of interest.
#' - If \code{do_pca} is TRUE, other results from principal component analysis.
#'
#' @usage
#' master_matrix(region, variables, do_pca = FALSE, center = TRUE, scale = FALSE,
#'               variables_in_matrix = NULL)
#'
#' @export
#' @importFrom raster mask crop rasterToPoints
#' @importFrom stats prcomp
#'
#' @examples
#' # Data
#' data("mx", package = "biosurvey")
#' variables <- raster::stack(system.file("extdata/variables.tif",
#'                                        package = "biosurvey"))
#'
#' # Create master matrix object
#' m_matrix <- master_matrix(region = mx, variables = variables, do_pca = TRUE,
#'                           center = TRUE, scale = TRUE)


master_matrix <- function(region, variables, do_pca = FALSE, center = TRUE,
                          scale = FALSE, variables_in_matrix = NULL) {
  # Initial tests
  if (missing(region)) {
    stop("Argument 'region' must be defined")
  }

  if (class(region)[1] != "SpatialPolygonsDataFrame") {
    stop("'region' must be of class SpatialPolygonsDataFrame")
  }

  if (missing(variables)) {
    stop("Argument 'variables' must be defined")
  }

  if (!class(variables)[1] %in% c("RasterStack", "RasterBrick")) {
    stop("'variables' must be of class RasterStack or RasterBrick")
  }

  # Mask variables to polygons (region of interest)
  variables <- raster::mask(raster::crop(variables, region), region)

  # Raster to matrix
  variables <- raster::rasterToPoints(variables)

  # Selecting variables for the matrix
  if (is.null(variables_in_matrix)) {
    vars <- variables[, -(1:2)]
  } else {
    vars <- variables[, variables_in_matrix]
  }

  # If do_pca is TRUE, do PCA
  if (do_pca == TRUE) {
    pca <- stats::prcomp(variables[, -(1:2)], center = center, scale. = scale)

    # Create matrix
    master_m <- data.frame(Longitude = variables[, 1], Latitude = variables[, 2],
                           vars, pca$x[, 1:2])

    # Return results
    return(list(master_matrix = master_m, polygon = region, PCA_results = pca))

  } else {
    # Create matrix
    master_m <- data.frame(Longitude = variables[, 1], Latitude = variables[, 2],
                           vars)

    # Return results
    return(list(master_matrix = master_m, polygon = region))
  }
}
