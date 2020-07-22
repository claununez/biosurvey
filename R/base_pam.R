#' Prepare a presence-absence matrix (PAM)
#'
#' @description Prepare a presence-absence matrix (PAM) in which all initial points
#' of interest (rows) will have a representation of the species present in such
#' areas (columns). Initial points of interest will be represented by an ID, and
#' longitude and latitude coordinates.
#'
#' @param data species geographic ranges to be used to create a presence-absence
#' matrix (PAM). This argument can be: data.frame, RasterStack, RasterBrick,
#' list, SpatialPolygonsDataFrame, SpatialPointsDataFrame, or character. See
#' details for description of characteristics of each option.
#' @param format (character) if \code{data} is a character, available formats are:
#' "shp", "gpkg", "GTiff", and "ascii".
#' @param master_matrix object of class "master_matrix" or "master_selection". See
#' details.
#' @param cell_size (numeric) resolution for grid (single number or vector of two
#' numbers) in kilometers (km).
#' @param complete_cover (logical) whether or not to include cells of grid
#' partially overlapped with the geographic region of interest contained in
#' \code{master_matrix}. Default = TRUE.
#'
#' @details
#' Objects of class "master_matrix" or "master_selection" could be obtained from
#' functions \code{\link{master_matrix}}, \code{\link{random_selection}},
#' \code{\link{uniformG_selection}}, \code{\link{uniformE_selection}}, or
#' \code{uniformEG_selection}.
#'
#' Description of objects to be used as \code{data}:
#' - data.frame.- a table containing  three columns. Columns must be in the
#' following order: Longitude, Latitude, Species.
#' - RasterStack or RasterBrick.- Each layer must be named as the species which
#' range it represents, and values in each layer must be 1 (presence) and 0
#' (absence).
#' - list.-  a list of RasterLayers that cannot be stacked because of extent or
#' resolution differences. Each element of the list must be named as the species
#' which range it represents, and values in each RasterLayer must be 1 (presence)
#' and 0 (absence).
#' - SpatialPolygonsDataFrame.- object representing species' geographic ranges.
#' The data frame associated with the object must contain a column
#' named "Species" to distinguish among features representing each species range.
#' - SpatialPointsDataFrame.- object of spatial points where each record of a
#' species must be a point. The associated data.frame must contain the following
#' columns (in that order): Longitude, Latitude, Species.
#' - character.- name of directory containing raster, shapefiles, or geopackage
#' files representing species geographic ranges. Each file must be named as the
#' species that it represents. All files must be in an only format. If files are
#' raster, values in each layer must be 1 (presence) and 0 (absence).
#'
#' @return
#' A presence-absence matrix (PAM) of class base_PAM for the region of interest
#' associated with a SpatialPolygonsDataFrame, as in a grid of \code{cell_size}
#' resolution. Each grid cell is related to a specific ID and longitude and
#' latitude coordinates. Presence (1) and absence (0) values for each species
#' in every cell of the PAM are included as apart of the data frame of the
#' SpatialPolygonsDataFrame. PAM statistics is returned as NULL, but other
#' functions can be used to calculate such statistics.
#'
#' @usage
#' base_PAM(data, format = NULL, master_matrix, cell_size,
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
#' # Create base_PAM
#' b_pam <- base_PAM(data = species_data, master_matrix = m_matrix, cell_size = 100)
#' sp::plot(b_pam$PAM)
#' summary(b_pam$PAM@data[, 1:6])


base_PAM <- function(data, format = NULL, master_matrix, cell_size,
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

  # Where to prepare spatial PAM
  where <- ifelse(!is.null(master_matrix$mask), "mask", "region")

  # Create geographic grid
  grid_r_pol <- grid_from_region(master_matrix[[where]], cell_size, complete_cover)

  # Prepare SpaptialPoints from different objects
  if (!is.data.frame(data)) {
    ## from raster objects
    if (clsdata %in% c("RasterStack", "RasterBrick")) {
      data <- stack_2data(species_layers = data)
    }

    ## from a list
    if (clsdata == "list") {
      data <- rlist_2data(raster_list = data)
    }

    ## from files stored in a directory
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

  # Assign ID to points depending of type of data if needed
  ## direct step from SPDF or character data argument
  if (clsdata %in% c("SpatialPolygonsDataFrame", "character")) {
    if (clsdata == "SpatialPolygonsDataFrame") {
      sp_points <- spdf_2data(spdf_object = data, spdf_grid = grid_r_pol)
    }
    if (clsdata == "character" & clsdata %in% c("shp", "gpkg")) {
      sp_points <- files_2data(path = data, format, spdf_grid = grid_r_pol)
    }
  }

  ## merging with ID if any other object is defined in data
  if (class(data)[1] %in% c("data.frame", "SpatialPointsDataFrame")) {
    if (class(data)[1] %in% "SpatialPointsDataFrame") {
      sp_points <- data
    }
    sp_points <- data.frame(ID = sp::over(methods::as(sp_points, "SpatialPoints"),
                                          grid_r_pol[, "ID"]),
                            Species = sp_points@data[, 3])
  }

  # PAM from points
  sp_points <- PAM_from_table(sp_points, ID_column = "ID",
                              species_column = "Species")

  # Complete PAM
  grid_r_pol@data <- merge(grid_r_pol@data, sp_points, by = "ID", all.x = TRUE)
  grid_r_pol@data[is.na(grid_r_pol@data)] <- 0

  # Returning results
  return(structure(list(PAM = grid_r_pol, PAM_indices = NULL), class = "base_PAM"))
}


