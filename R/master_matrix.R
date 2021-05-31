#' Prepare a base object to perform further analyses
#'
#' @description prepare an S3 object that will serve as the base to perform all
#' further analyses. This object will contain geographic and environmental
#' information that will be used to characterize the region of interest.
#'
#' @param region SpatialPolygons* of the region of interest; for instance,
#' a country, another type of administrative are, or a protected area.
#' @param variables RasterStack or RasterBrick of environmental variables.
#' @param mask (optional) SpatialPolygons* object to mask \code{variables} and
#' reduce \code{region} to an area that is more relevant for analysis (e.g.,
#' only areas with natural vegetation cover). Default = NULL.
#' @param preselected_sites data.frame containing sites that must be included
#' in posterior selections of sites for the survey system. Columns must be:
#' "Sites", "Longitude", "Latitude", in that order.
#' @param do_pca (logical) whether or not to perform a principal component
#' analysis. Default = FALSE.
#' @param center (logical) whether or not to center variables. Argument to be
#' passed to the function \code{\link[stats]{prcomp}}. Default = TRUE.
#' @param scale (logical) whether or not to scale the variables. Recommended
#' when variables are in different units. Argument to be passed to the function
#' \code{\link[stats]{prcomp}}. Default = FALSE.
#' @param variables_in_matrix (character) name of variables to include in
#' matrix. If NULL (the default) all variables will be included.
#' @param verbose (logical) whether or not to print messages about the process.
#' Default = TRUE.
#'
#' @return
#' An S3 object of class \code{\link{master_matrix}} containing the following
#' elements:
#' - data_matrix: a date.frame with information about geographic location of
#' raster cells, initial environmental data, and if \code{do_pca} is TRUE,
#' the first two principal components derived from original data.
#' - region: a SpatialPolygons* representing the region of interest.
#' - mask: SpatialPolygons* object used. NULL if \code{mask} was not defined.
#' - preselected_sites: sites defined by used. NULL if \code{preselected_sites}
#' was not defined.
#' - raster_base: a SpatialPolygonsDataFrame representing the grid of the
#' rasters used, which will be used for plotting purposes.
#' - PCA_results: if \code{do_pca} is TRUE, other results from principal
#' component analysis. If FALSE, PCA_results element of the object is NULL.
#'
#' @details
#' This function helps in preparing all data as needed for further analyses
#' aiming to define a survey sampling system considering geographic and
#' environmental spaces in the \code{region} of interest.
#'
#' If \code{mask} is defined all analyses will be restricted to such an area.
#' If \code{mask} is not fully contained by \code{region}, the mask used for
#' reducing \code{variables}, and returned as part of the S3 object
#' (master_matrix) is the intersection between them.
#'
#' If \code{preselected_sites} is defined, environmental values and, if
#' \code{do_pca} = TRUE, principal components are added to such records. These
#' records and their characteristics will be considered in further analyses.
#'
#' @usage
#' prepare_master_matrix(region, variables, mask = NULL,
#'                       preselected_sites = NULL, do_pca = FALSE,
#'                       center = TRUE, scale = FALSE,
#'                       variables_in_matrix = NULL, verbose = TRUE)
#'
#' @export
#' @importFrom raster mask crop rasterToPoints intersect
#' @importFrom stats prcomp predict
#'
#' @examples
#' # Data
#' data("mx", package = "biosurvey")
#' variables <- raster::stack(system.file("extdata/variables.tif",
#'                                        package = "biosurvey"))
#'
#' # Create master matrix object
#' m_matrix <- prepare_master_matrix(region = mx, variables = variables,
#'                                   do_pca = TRUE, center = TRUE, scale = TRUE)


prepare_master_matrix <- function(region, variables, mask = NULL,
                                  preselected_sites = NULL, do_pca = FALSE,
                                  center = TRUE, scale = FALSE,
                                  variables_in_matrix = NULL, verbose = TRUE) {
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
    #mask <- raster::intersect(region, mask)
  }

  # Mask variables to polygons (region of interest)
  if (verbose == TRUE) {
    message("Processing raster layers, please wait...")
  }
  if (!is.null(mask)) {
    variables <- raster::mask(raster::crop(variables, mask), mask)
  } else {
    variables <- raster::mask(raster::crop(variables, region), region)
  }

  # Adding user defined points
  if (!is.null(preselected_sites)) {
    if (class(preselected_sites)[1] != "data.frame") {
      stop("Argument 'preselected_sites' must be of class data.frame.")
    }

    preselected_sites <- data.frame(preselected_sites,
                                    raster::extract(variables,
                                                    preselected_sites[, 2:3]))
    colnames(preselected_sites)[2:3] <- c("Longitude", "Latitude")
  }

  # Base raster
  b_raster <- variables[[1]]
  names(b_raster) <- "base"
  b_raster[!is.na(b_raster[])] <- 1
  b_raster <- as(b_raster,"SpatialPolygonsDataFrame")

  # Raster to matrix
  variables <- raster::rasterToPoints(variables)
  colnames(variables)[1:2] <- c("Longitude", "Latitude")

  # Selecting variables for the matrix
  if (!is.null(variables_in_matrix)) {
    variables <- variables[, c("Longitude", "Latitude", variables_in_matrix)]
    if (!is.null(preselected_sites)) {
      prna <- colnames(preselected_sites)
      preselected_sites <- preselected_sites[, c(prna[1:3], variables_in_matrix)]
    }
  }

  # If do_pca is TRUE, do PCA
  if (do_pca == TRUE) {
    if (verbose == TRUE) {
      message("Performing PCA analysis")
    }
    pca <- stats::prcomp(variables[, -(1:2)], center = center, scale. = scale)

    # Create matrix
    master_m <- data.frame(variables, pca$x[, 1:2])

    # Predict in predefined sites
    if (!is.null(preselected_sites)) {
      sspca <- stats::predict(pca, preselected_sites[, -(1:3)])
      preselected_sites <- data.frame(preselected_sites, sspca[, 1:2])
    }

    # Return results
    return(new_master_matrix(data_matrix = master_m, region = region,
                             mask = mask, raster_base = b_raster,
                             preselected_sites = preselected_sites,
                             PCA_results = pca))

  } else {
    # Create matrix
    master_m <- data.frame(variables)

    # Return results
    return(new_master_matrix(data_matrix = master_m, region = region,
                             mask = mask, raster_base = b_raster,
                             preselected_sites = preselected_sites,
                             PCA_results = NULL))
  }
}
