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
