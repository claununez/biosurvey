#' Presence-absence matrix (PAM) linked to a spatial grid
#'
#' @description Prepares a presence-absence matrix (PAM) in which all sites
#' of interest (rows) will have a value for presence or absence of a species
#' of interest (columns). Initial points of interest will be represented by an
#' ID, and longitude and latitude coordinates. The PAM will be linked to a
#' spatial grid.
#'
#' @param data species geographic ranges to be used to create a presence-absence
#' matrix (PAM). This argument can be: character, data.frame, SpatRaster,
#' list, or SpatVector. See details for a description of the characteristics of
#' data for each option.
#' @param format (character) if \code{data} is of class character, available
#' options are: "shp", "gpkg", "GTiff", and "ascii".
#' @param master_matrix object of class "master_matrix" or "master_selection".
#' See details. Default = NULL. Either this argument or \code{region} must be
#' defined.
#' @param region SpatVector of the region of interest; for instance,
#' a country, another type of administrative are, or a protected area.
#' Default = NULL. Either this argument or \code{master_matrix} must be
#' defined. Ignored if \code{master_matrix} is defined.
#' @param cell_size (numeric) resolution for grid (single number or vector of
#' two numbers) in kilometers (km).
#' @param complete_cover (logical) whether or not to include cells of grid
#' partially overlapped with the geographic region of interest contained in
#' \code{master_matrix}. Default = TRUE.
#' @param clip_grid (logical) whether to clip the spatial grid using the region
#' of interest. Clipping improves visualization but depending on how complex
#' the region of interest is it could take time to perform this task.
#' @param indices (character) code for indices to be calculated. Basic indices
#' are calculated all the time, other indices need to be specified. Options are:
#' "all", "basic, "AB", "BW", "BL", "SCSC", "SCSR", "DF", "DivF", "CC", "WRN",
#' "SRC", "CMSC", and "CMSR". Default = "basic". See details.
#' @param parallel (logical) whether to perform analyses in parallel.
#' Default = FALSE. Not used if data is of class data.frame, or
#' SpatRaster.
#' @param n_cores (numeric) number of cores to be used when \code{parallel} =
#' TRUE. The default, NULL, uses available cores - 1.
#' @param verbose (logical) whether or not to print messages about the process.
#' Default = TRUE.
#'
#' @details
#' Objects of class "master_matrix" or "master_selection" can be obtained from
#' functions \code{\link{prepare_master_matrix}},
#' \code{\link{random_selection}}, \code{\link{uniformG_selection}},
#' \code{\link{uniformE_selection}}, or \code{uniformEG_selection}. The element
#' region or mask if this last is not NULL is used to prepare the spatial grid.
#'
#' Geographic projection of objects or coordinates involved must be WGS84
#' (EPSG:4326).
#'
#' Description of objects to be used as \code{data}:
#' - character.- name of directory containing raster, shapefiles or geopackage
#' files representing species geographic ranges. Each file must be named
#' as the species that it represents. All files must be in the same format. If
#' files are in raster format, "GTiff" and "ascii" are acceptable extensions;
#' values in each layer must be 1 (presence) and 0 (absence).
#' - data.frame.- a table containing  three columns. Columns must be in the
#' following order: Longitude, Latitude, Species.
#' - SpatRaster.- Each layer must be named as the species which
#' range it represents, and values in each layer must be 1 (presence) and 0
#' (absence).
#' - list.-  a list of SpatRaster objects that cannot be stacked because of
#' extent or resolution differences. Each element of the list must be named as
#' the species which range it represents, and values in each SpatRaster must be
#' 1 (presence, suitable) and 0 (absence, unsuitable).
#' - SpatVector.- object representing species' geographic ranges.
#' The data.frame associated with the object must contain a column named
#' "Species" to distinguish among features representing each species range.
#' - SpatVector.- object of spatial points where each record of a
#' species must be a point. The associated data.frame must contain the
#' following columns (in that order): Longitude, Latitude, Species.
#'
#' A list of codes and indices that can be calculated is described below. For
#' further details on the way calculations are performed and the meaning of the
#' indices see Soberon and Cavner (2015)
#' \doi{https://doi.org/10.17161/bi.v10i0.4801}.
#'
#' |Code  |Index                                    |Calculation                     |
#' |:-----|----------------------------------------:|-------------------------------:|
#' |RI    |Richness                                 |Basic                           |
#' |RA    |Range                                    |Basic                           |
#' |RIN   |Richness normalized                      |Basic                           |
#' |RAN   |Range normalized                         |Basic                           |
#' |AB    |Additive Beta                            |Needs to be defined             |
#' |BW    |Beta Whittaker                           |Needs to be defined             |
#' |BL    |Beta Legendre                            |Needs to be defined and DF      |
#' |SCSC  |Schluter covariance sites-composition    |Needs to be defined and CMSC    |
#' |SCSR  |Schluter covariance species-ranges       |Needs to be defined and CMSR    |
#' |DF    |Dispersion field                         |Needs to be defined             |
#' |DivF  |Diversity field                          |Needs to be defined             |
#' |SCC   |Shared community composition             |Needs to be defined             |
#' |WRN   |Wright-Reeves nestedness                 |Needs to be defined, BW, and DF |
#' |SRC   |Stone-Roberts C-score                     |Needs to be defined and DF      |
#' |CMSC  |Covariance matrix sites-composition      |Needs to be defined, DF, and BW |
#' |CMSR  |Covariance matrix species-ranges         |Needs to be defined, SCC, and BW|
#' |MCC   |Mean composition covariance              |Calculated with CMSC            |
#' |MRC   |Mean range covariance                    |Calculated with CMSR            |
#'
#' @seealso \code{\link{PAM_indices}}
#'
#' @return
#' A presence-absence matrix (PAM) of class \code{\link{base_PAM}} for the
#' region of interest associated with a SpatVector, as in a grid
#' of \code{cell_size} resolution. Each grid cell is related to a specific ID
#' and longitude and latitude coordinates. Presence (1) and absence (0) values
#' for each species in every cell of the PAM are included as apart of the
#' data.frame of the SpatVector. PAM indices is returned with the
#' basic indices of biodiversity as default, but can be changed using the
#' argument \code{indices}.
#'
#' @usage
#' prepare_base_PAM(data, format = NULL, master_matrix = NULL, region = NULL,
#'                  cell_size, complete_cover = TRUE, clip_grid = FALSE,
#'                  indices = "basic", parallel = FALSE, n_cores = NULL,
#'                  verbose = TRUE)
#'
#' @export
#' @importFrom terra vect extract crs values intersect
#'
#' @examples
#' # Data
#' m_matrix <- read_master(system.file("extdata/m_matrix.rds",
#'                                     package = "biosurvey"))
#' species_data <- terra::vect(system.file("extdata/species_data.gpkg",
#'                                         package = "biosurvey"))
#'
#' # Create base_PAM
#' b_pam <- prepare_base_PAM(data = species_data, master_matrix = m_matrix,
#'                           cell_size = 100)
#' terra::plot(b_pam$PAM)
#' summary(b_pam)

prepare_base_PAM <- function(data, format = NULL, master_matrix = NULL,
                             region = NULL, cell_size, complete_cover = TRUE,
                             clip_grid = FALSE, indices = "basic",
                             parallel = FALSE, n_cores = NULL, verbose = TRUE) {
  # Initial tests
  if (missing(data)) {
    stop("Argument 'data' must be defined")
  }
  if (is.null(master_matrix) & is.null(region)) {
    stop("One of the arguments 'master_matrix' or 'region' must be defined.")
  }

  clsdata <- class(data)[1]

  if (!clsdata %in% c("SpatRaster", "data.frame", "list",
                      "SpatVector", "character")) {
    stop("Argument 'data' is not valid, check function's help.")
  }

  if (clsdata == "character") {
    if (dir.exists(data) == FALSE) {
      stop("Directory defined in 'data' not found")
    }
    if (is.null(format)) {
      stop("Argument 'format' must be defined if class of 'data' is character.")
    }
  }

  all_in <- c("all", "basic", "AB", "BW", "BL", "SCSC", "SCSR", "DF", "SCC",
              "WRN", "SRC", "CMSC", "CMSR")
  if (any(!indices %in% all_in)) {
    stop("One or more elements defined in 'indices' is not valid, check function's help.")
  }

  # Create geographic grid
  if (verbose == TRUE) {
    message("Preparing spatial grid")
  }

  if (!is.null(master_matrix)) {
    # Where to prepare spatial PAM
    where <- ifelse(!is.null(master_matrix$mask), "mask", "region")

    grid_r_pol <- grid_from_region(master_matrix[[where]], cell_size,
                                   complete_cover)
  } else {
    grid_r_pol <- grid_from_region(region, cell_size, complete_cover)
  }

  # Prepare SpatVector from different objects
  if (verbose == TRUE) {
    message("Preprocessing 'data'")
  }
  if (clsdata != "data.frame") {
    ## From raster objects
    if (clsdata %in% "SpatRaster") {
      data <- stack_2data(species_layers = data)
    }

    ## From a list
    if (clsdata == "list") {
      data <- rlist_2data(raster_list = data)
    }

    ## From files stored in a directory
    if (clsdata == "character") {
      if (!format %in% c("shp", "gpkg")) {
        data <- files_2data(path = data, format = format, parallel = parallel,
                            n_cores = n_cores)
      }
    }
  } else {
    sp_points <- terra::vect(data, geom = colnames(data)[1:2],
                             crs = terra::crs("EPSG:4326"))
  }

  # SpatialVector from data if needed
  if (!clsdata %in% "SpatVector") {
    if (clsdata == "character") {
      if (!format %in% c("shp", "gpkg")) {
        sp_points <- terra::vect(data, geom = colnames(data)[1:2],
                                 crs = terra::crs("EPSG:4326"))
      }
    } else if (clsdata %in% c("SpatRaster", "list")) {
      sp_points <- terra::vect(data, geom = colnames(data)[1:2],
                               crs = terra::crs("EPSG:4326"))
    }
  }

  # Assign ID to points depending of type of data if needed
  ## Direct step from SpatVector or character data argument
  if (clsdata %in% c("SpatVector", "character")) {
    if (clsdata == "SpatVector") {
      sp_points <- spdf_2data(spdf_object = data, spdf_grid = grid_r_pol,
                              parallel = parallel, n_cores = n_cores)
    }
    if (clsdata == "character") {
      if (format %in% c("shp", "gpkg")) {
        sp_points <- files_2data(path = data, format, spdf_grid = grid_r_pol,
                                 parallel = parallel, n_cores = n_cores)
      }
    }
  }

  ## Merging with ID if any other object is defined in data
  if (class(data)[1] == "data.frame") {
    sp_points <- data.frame(ID = terra::extract(grid_r_pol, sp_points)$ID,
                            Species = sp_points[[1]])
    sp_points <- unique(sp_points)
  }

  # PAM from points
  if (verbose == TRUE) {
    message("\nPreparing PAM from information")
  }
  sp_points <- PAM_from_table(sp_points, ID_column = "ID",
                              species_column = "Species")

  # Complete PAM
  terra::values(grid_r_pol) <- merge(grid_r_pol, sp_points, by = "ID",
                                     all.x = TRUE)
  terra::values(grid_r_pol)[is.na(terra::values(grid_r_pol))] <- 0
  coltk <- names(grid_r_pol)

  # Clipping if needed
  if (clip_grid == TRUE) {
    if (verbose == TRUE) {
      message("Clipping PAM to region of interest")
    }

    if (!is.null(master_matrix)) {
      grid_r_pol <- terra::intersect(grid_r_pol, master_matrix[[where]])
    } else {
      grid_r_pol <- terra::intersect(grid_r_pol, region)
    }
    coltk <- make.names(coltk)
    coltk <- intersect(coltk, names(grid_r_pol))
    terra::values(grid_r_pol) <- terra::values(grid_r_pol)[, coltk]
  }

  # Preparing and returning results
  if (verbose == TRUE) {
    message("Calculating PAM indices")
  }
  bPAM <- new_base_PAM(PAM = grid_r_pol, PAM_indices = NULL)

  bPAM <- PAM_indices(PAM = bPAM, indices = indices, exclude_column = NULL)

  return(bPAM)
}
