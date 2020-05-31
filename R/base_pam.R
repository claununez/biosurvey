#' Prepare a presence-absence matrix (PAM)
#'
#' @description Prepare a presence-absence matrix (PAM) in which all initial points
#' of interest (rows) will have a representation of the species present in such
#' areas (columns). Initial points of interest will be represented by an ID, and
#' longitude and latitude coordinates.
#'
#' @param data species geographic ranges to be used to create a presence-absence
#' matrix (PAM). This argument can be: RasterStack, RasterBrick, data.frame,
#' list, SpatialPolygonsDataFrame, SpatialPointsDataFrame, or character. See
#' details for description of characteristics of each option.
#' @param format (character) if \code{data} is a character, available formats are:
#' "shp", "gpkg", "GTiff", and "ascii".
#' @param master_matrix object derived from function \code{\link{master_matrix}}.
#' Optionally, if master_matrix is not necessary, a list containing an object of
#' class SpatialPolygonsDataFrame, representinng the region of interest, can be
#' used. The name of this element in the list must be "polygon". For instance:
#' \code{my_list <- list(polygon = YOUR_SpatialPolygonsDataFrame)}.
#' @param cell_size (numeric) resolution for grid (single number or vector of two
#' numbers) in kilometers (km).
#' @param complete_cover (logical) whether or not to include cells of grid
#' partially overlapped with the geographic region of interest contined in
#' \code{master_matrix}. Default = TRUE.
#'
#' @return
#' A presence-absence matrix (PAM) for the region of interest associated with a
#' SpatialPolygonsDataFrame, as in a grid of \code{cell_size} resolution. Each
#' grid cell is related to a specific ID and longitude and latitude coordinates.
#' Presence (1) and absence (0) values for each species in every cell of the PAM
#' are included as apart of the data frame of the SpatialPolygonsDataFrame.
#'
#' @usage
#' base_pam(data, format = NULL, master_matrix, cell_size,
#'          complete_cover = TRUE)
#'
#' @export
#' @importFrom sp SpatialPointsDataFrame over
#' @importFrom methods as
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#' data("species_data", package = "biosurvey")
#'
#' # Create base_pam
#' b_pam <- base_pam(data = species_data, master_matrix = m_matrix, cell_size = 100)
#'
#' sp::plot(b_pam)
#' summary(b_pam@data[, 1:6])

base_pam <- function(data, format = NULL, master_matrix, cell_size,
                     complete_cover = TRUE) {
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
  grid_r_pol <- grid_from_region(master_matrix$polygon, cell_size, complete_cover)

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
    if (class(data)[1] %in% "SpatialPointsDataFrame") {
      sp_points <- data
    }
    sp_points <- data.frame(ID = sp::over(methods::as(sp_points, "SpatialPoints"),
                                          grid_r_pol[, "ID"]),
                            Species = sp_points@data[, 3])
  }

  # PAM from points
  sp_points <- pam_from_table(sp_points, ID_column = "ID",
                              species_column = "Species")

  # Complete PAM
  grid_r_pol@data <- merge(grid_r_pol@data, sp_points, by = "ID", all.x = TRUE)
  grid_r_pol@data[is.na(grid_r_pol@data)] <- 0

  return(grid_r_pol)
}


