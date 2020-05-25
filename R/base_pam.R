#' Prepare a presence-absence matrix (PAM)
#'
#' @description Prepare a presence-absence matrix (PAM) in which all initial points
#' of interest (rows) will have a representation of the species present in such
#' areas (columns). Initial points of interest are represented by an ID, and
#' longitude and latitude coordinates.
#'
base_pam <- function(data, format = NULL, master_matrix, cell_size) {
  # Initial tests
  clsdata <- class(data)[1]

  if (!clsdata %in% c("RasterStack", "RasterBrick", "data.frame", "list",
                     "SpatialPolygonsDataFrame", "SpatialPointsDataFrame",
                     "character")) {
    stop("Argument 'data' is not valid, check function's help.")
  }

  if (clsdata == "character") {
    if (is.null(format)) {
      stop("Argument 'format' must be defined if class of 'data' is character")
    }
  }

  # Create geographyc grid
  grid_r_pol <- grid_from_region(master_matrix$polygon, cell_size)

  # Prepare SpaptialPoints from different objects
  if (!is.data.frame(data)) {
    if (clsdata %in% c("RasterStack", "RasterBrick")) {
      data <- stack_2data(species_layers = data)
    }

    if (clsdata == "list") {
      data <- rlist_2data(raster_list = data)
    }

    if (clsdata == "character" & !clsdata %in% c("shp", "gpkg")) {
      data <- files_2data(path = data, format)
    }
  }

  # SpatialPointsDataFrame from data if needed
  if (!clsdata %in% c("SpatialPointsDataFrame", "SpatialPolygonsDataFrame")) {
    if (clsdata == "character" & !clsdata %in% c("shp", "gpkg")) {
      sp_points <- sp::SpatialPointsDataFrame(data[, 1:2], data = data)
    }
  }

  # Asign ID to points
  if (clsdata %in% c("SpatialPolygonsDataFrame", "character")) {
    if (clsdata == "SpatialPolygonsDataFrame") {
      sp_points <- spdf_2data(spdf_object = data, spdf_grid = grid_r_pol)
    }
    if (clsdata == "character" & clsdata %in% c("shp", "gpkg")) {
      sp_points <- files_2data(path = data, format, spdf_grid = grid_r_pol)
    }
  }

  if (class(data)[1] %in% c("data.frame", "SpatialPointsDataFrame")) {
    sp_points <- data.frame(ID = sp::over(as(sp_points, "SpatialPoints"),
                                          grid_r_pol[, "ID"]),
                            Species = sp_points@data[, 3])
  }

  # PAM from points
  sp_points <- pam_from_table(sp_points, id_column = "ID",
                              species_column = "Species")

  # Complete PAM
  grid_r_pol@data <- merge(grid_r_pol@data, sp_points, by = "ID", all.x = TRUE)
  grid_r_pol@data[is.na(grid_r_pol@data)] <- 0

  return(grid_r_pol)
}


