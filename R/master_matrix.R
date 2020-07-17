#' Prepare a base object to perform further analyses
#'
#' @description Prepare an S3 object that will serve as the base to perform all
#' further analyses. This object will contain geographic and environmental
#' information that will be used to characterize the region of interest.
#'
#' @param region SpatialPolygons* of the region of interest; for instance,
#' a country, another type of administrative are, or a protected area.
#' @param variables RasterStack or RasterBrick of environmental variables.
#' @param mask (optional) SpatialPolygons* object to mask \code{variables} and
#' reduce \code{region} to an area that is more relevant for analysis (e.g., only
#' areas with natural vegetation cover). Default = NULL.
#' @param do_pca (logical) whether or not to perform a Principal Component Analysis.
#' Default = FALSE.
#' @param center (logical) whether or not to center variables. Argument to be passed
#' to the function \code{\link[stats]{prcomp}}. Default = TRUE.
#' @param scale (logical) whether or not to scale the variables. Recommended when
#' variables are in different units. Argument to be passed to the function
#' \code{\link[stats]{prcomp}}. Default = FALSE.
#' @param variables_in_matrix (character) name of variables to include in matrix.
#' If NULL (the default) all variables will be included.
#'
#' @return
#' An S3 object of class master_matrix containing the following elements:
#' - master_matrix: a date.frame with information about geographic location of
#' raster cells, initial environmental data, and if \code{do_pca} is TRUE,
#' the first two principal components derived from original data.
#' - region: a SpatialPolygons* representing the region of interest.
#' - raster_base: a raster layer for the region of interest with a single value,
#' to be used for plotting purposes.
#' - mask: SpatialPolygons* object used. Present only if \code{mask} was defined.
#' - PCA_results: if \code{do_pca} is TRUE, other results from principal
#' component analysis. If FALSE, PCA_results element of the object is NULL.
#'
#' @details
#' This function helps in preparing all data as needed for further analyses
#' aiming to define a survey sampling system considering geographic and
#' environmental spaces in the \code{region} of interest.
#'
#' If \code{mask} is defined all analyses will be restricted to such an area. If
#' \code{mask} is not fully contained by \code{region}, the mask used for reducing
#' \code{variables}, and returned as part of the S3 object (master_matrix) is the
#' intersection between them.
#'
#' @usage
#' master_matrix(region, variables, mask = NULL, do_pca = FALSE, center = TRUE,
#'               scale = FALSE, variables_in_matrix = NULL)
#'
#' @export
#' @importFrom raster mask crop rasterToPoints intersect
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


master_matrix <- function(region, variables, mask = NULL, do_pca = FALSE,
                          center = TRUE, scale = FALSE, variables_in_matrix = NULL) {
  # Initial tests
  if (missing(region)) {
    stop("Argument 'region' must be defined")
  }
  if (!class(region)[1] %in% c("SpatialPolygons", "SpatialPolygonsDataFrame")) {
    stop("'region' must be of class SpatialPolygons*")
  }
  if (missing(variables)) {
    stop("Argument 'variables' must be defined")
  }
  if (!class(variables)[1] %in% c("RasterStack", "RasterBrick")) {
    stop("'variables' must be of class RasterStack or RasterBrick")
  }
  if (!is.null(mask)) {
    if (!class(mask)[1] %in% c("SpatialPolygons", "SpatialPolygonsDataFrame")) {
      stop("'mask' must be of class SpatialPolygons*")
    }

    # Intersection of region and mask to avoid complications
    mask <- raster::intersect(region, mask)
  }

  # Mask variables to polygons (region of interest)
  if (!is.null(mask)) {
    variables <- raster::mask(raster::crop(variables, mask), mask)
  } else {
    variables <- raster::mask(raster::crop(variables, region), region)
  }

  # Base raster
  b_raster <- variables[[1]]
  names(b_raster) <- "base"
  b_raster[!is.na(b_raster[])] <- 1

  # Raster to matrix
  variables <- raster::rasterToPoints(variables)
  colnames(variables)[1:2] <- c("Longitude", "Latitude")

  # Selecting variables for the matrix
  if (!is.null(variables_in_matrix)) {
    variables <- variables[, c("Longitude", "Latitude", variables_in_matrix)]
  }

  # If do_pca is TRUE, do PCA
  if (do_pca == TRUE) {
    pca <- stats::prcomp(variables[, -(1:2)], center = center, scale. = scale)

    # Create matrix
    master_m <- data.frame(variables, pca$x[, 1:2])

    # Return results
    return(structure(list(master_matrix = master_m, polygon = region, mask = mask,
                          raster_base = b_raster, PCA_results = pca),
                     class = "master_matrix"))

  } else {
    # Create matrix
    master_m <- data.frame(variables)

    # Return results
    return(structure(list(master_matrix = master_m, region = region, mask = mask,
                          raster_base = b_raster, PCA_results = NULL),
                     class = "master_matrix"))
  }
}
