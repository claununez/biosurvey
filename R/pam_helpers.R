#' Creates grid for a given geographic region
#'
#' @description Divides the region of interest in a grid of a specific cell size.
#'
#' @param region SpatVector object of a polygon for the region of interest.
#' Object projection must be World Geodetic System (WGS84).
#' @param cell_size (numeric) resolution for grid (single number or vector of
#' two numbers) in kilometers (km).
#' @param complete_cover (logical) whether or not to include cells of grid
#' partially overlapped with region. Default = TRUE.
#'
#' @return
#' Grid SpatVector for the region of interest. Each grid cell
#' is related to a specific ID and longitude and latitude coordinates.
#'
#' @usage
#' grid_from_region(region, cell_size, complete_cover = TRUE)
#'
#' @export
#' @importFrom terra mask rasterize crs ext project geom rast as.data.frame
#'
#' @examples
#' # Data
#' mx <- terra::vect(system.file("extdata/mx.gpkg", package = "biosurvey"))
#'
#' # Create grid from polygon
#' grid_reg <- grid_from_region(region = mx, cell_size = 100)
#'
#' terra::plot(grid_reg)

grid_from_region <- function(region, cell_size, complete_cover = TRUE) {
  # Initial tests
  if (missing(region)) {
    stop("Argument 'region' must be defined")
  }
  if (class(region)[1] != "SpatVector") {
    stop("'region' must be of class 'SpatVector'")
  }
  if (missing(cell_size)) {
    stop("Argument 'cell_size' must be defined")
  } else {

    # Projecting region toLambert equeal area projection
    if (is.na(terra::crs(region))) {
      stop("'region' must be projected to WGS84 (EPSG:4326)")
    }
    WGS84 <- terra::crs("+init=epsg:4326")
    region <- terra::project(region, WGS84)
    cent <- terra::geom(terra::centroids(region))[, c("x", "y")]
    LAEA <- terra::crs(paste0("+proj=laea +lat_0=", cent[2], " +lon_0=",
                              cent[1], " +x_0=0 +y_0=0 +ellps=WGS84 ",
                              "+datum=WGS84 +units=m +no_defs"))
    region <- terra::project(region, LAEA)

    # Test if dimensions are valid
    dims <- terra::ext(region)
    xdim <- diff(dims[1:2])
    ydim <- diff(dims[3:4])
    if (length(cell_size) > 2) {
      stop("Argument 'cell_size' must be of length 1 or 2.")
    } else {
      if (length(cell_size) == 1) {
        cell_size <- rep(cell_size, 2)
      }
      if (cell_size[1] >= xdim & cell_size[2] >= ydim) {
        stop("'cell_size' must be smaller than at least one of the dimensions of 'region'")
      }
    }
  }

  # Creating a grid
  grid <- terra::rast(region, res = cell_size * 1000, vals = 1)

  # Extract grid with region
  grid_reg <- terra::mask(grid, region, touches = complete_cover)

  # Back to WGS84
  grid_reg <- terra::project(grid_reg, WGS84)

  # Grid for region of interest
  grid_r_pol <- terra::as.polygons(grid_reg, dissolve = FALSE)

  # Points for region of interest
  grid_r_pol <- cbind(grid_r_pol,
                      terra::as.data.frame(grid_reg, xy = TRUE,
                                           cells = TRUE)[, 1:3])
  grid_r_pol[, 1] <- NULL
  names(grid_r_pol) <- c("ID", "Longitude", "Latitude")

  return(grid_r_pol)
}



#' Creates a data.frame of species' references from SpatRaster
#'
#' @description Creates a data.frame of species' references that contains
#' longitude, latitude, and species name, using a SpatRaster representing
#' multiple species as input.
#'
#' @usage
#' stack_2data(species_layers)
#'
#' @param species_layers SpatRaster object. Each layer must be
#' named as the species that it represents, and values in each layer must be
#' 1 (presence) and 0 (absence).
#'
#' @return
#' A data.frame of species geographic records derived from values of presence
#' in each layer from the SpatRaster
#'
#' @export
#' @importFrom terra as.data.frame
#'
#' @examples
#' # Data
#' rsp <- terra::rast(system.file("extdata/sp_layers.tif",
#'                                 package = "biosurvey"))
#' names(rsp) <- paste0("Species_", 1:5)
#'
#' # Species data from SpatRaster
#' sp_data <- stack_2data(species_layers = rsp)
#' summary(sp_data)

stack_2data <- function(species_layers) {
  # Initial tests
  if (missing(species_layers)) {
    stop("Argument 'species_layers' must be defined")
  }
  if (class(species_layers)[1] != "SpatRaster") {
    stop("'species_layers' must be of class 'SpatRaster'")
  }

  # Stack to matrix
  sppm <- terra::as.data.frame(species_layers, xy = TRUE)
  spnames <- colnames(sppm)[-c(1, 2)]

  # Preparing data
  sps <- lapply(1:length(spnames), function(x) {
    cond <- sppm[, 2 + x] == 1
    data.frame(sppm[cond, 1], sppm[cond, 2], spnames[x])
  })

  sps <- do.call(rbind, sps)
  colnames(sps) <- c("Longitude", "Latitude", "Species")

  return(sps)
}



#' Creates a data.frame of species' references from SpatVector
#'
#' @description Creates a data.frame of species' references that contains
#' identifiers of position and species name, using a SpatVector representing
#' multiple species as input.
#'
#' @param spdf_object SpatVector representing species' geographic
#' distributions. The data.frame associated with the object must contain a
#' column named "Species" to distinguish among features.
#' @param spdf_grid SpatVector of geographic grid for the region of interest
#' (output of function \code{\link{grid_from_region}}).
#' @param parallel (logical) whether to perform analyses in parallel.
#' Default = FALSE.
#' @param n_cores (numeric) number of cores to be used when \code{parallel} =
#' TRUE. The default, NULL, uses available cores - 1.
#'
#' @return
#' A data.frame of species' found in distinct positions (defined with
#' identifiers); includes two columns: "ID" and "Species".
#'
#' @usage
#' spdf_2data(spdf_object, spdf_grid, parallel = FALSE, n_cores = NULL)
#'
#' @export
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom terra crs project
#' @importFrom stats na.omit
#'
#' @examples
#' # Data
#' species_data <- terra::vect(system.file("extdata/species_data.gpkg",
#'                                         package = "biosurvey"))
#' mx <- terra::vect(system.file("extdata/mx.gpkg", package = "biosurvey"))
#'
#' # GRID
#' grid_reg <- grid_from_region(region = mx, cell_size = 100)
#'
#' # Species data from polygons
#' sp_data <- spdf_2data(spdf_object = species_data, spdf_grid = grid_reg)
#' summary(sp_data)

spdf_2data <- function(spdf_object, spdf_grid, parallel = FALSE,
                       n_cores = NULL) {
  # Initial tests
  if (missing(spdf_object)) {
    stop("Argument 'spdf_object' must be defined")
  }
  if (missing(spdf_grid)) {
    stop("Argument 'spdf_grid' must be defined")
  }
  cond <- c(class(spdf_object)[1] != "SpatVector",
            class(spdf_grid)[1] != "SpatVector")
  if (any(cond)) {
    stop("'spdf_object' and 'spdf_grid' must be of class 'SpatVector'")
  }

  # Fixing projections
  if (terra::crs(spdf_object) != terra::crs(spdf_grid)) {
    spdf_object <- terra::project(spdf_object, terra::crs(spdf_grid))
  }

  # Names to be matched
  spnames <- as.character(spdf_object$Species)

  if (parallel == TRUE) {
    ## Preparing parallel running
    n_cores <- ifelse(is.null(n_cores), parallel::detectCores() - 1, n_cores)

    ## Progress combine (rbind) function
    pb <- utils::txtProgressBar(min = 1, max = length(spnames), style = 3)
    progress <- function(n) {
      utils::setTxtProgressBar(pb, n)
    }
    opts <- list(progress = progress)

    ## Make cluster
    cl <- snow::makeSOCKcluster(n_cores)
    doSNOW::registerDoSNOW(cl)

    ## wrap vectors
    spdf_grid <- terra::wrap(spdf_grid)
    spdf_object <- terra::wrap(spdf_object)

    ## Processing
    sps <- foreach::foreach(
      i = 1:length(spnames), .inorder = FALSE,
      .combine = "rbind", .options.snow = opts
    ) %dopar% {
      ID <- terra::unwrap(spdf_grid)[
        terra::unwrap(spdf_object)[spnames == spnames[i], ], ]

      return(na.omit(data.frame(ID = ID$ID, Species = spnames[i])))
    }

    snow::stopCluster(cl)
  } else {
    ## Progress bar
    pb <- utils::txtProgressBar(min = 1, max = length(spnames), style = 3)

    # Running in loop for all elements of list
    sps <- list()

    for (x in 1:length(spnames)) {
      Sys.sleep(0.1)
      utils::setTxtProgressBar(pb, x)

      # Preparing data
      ID <- spdf_grid[spdf_object[spnames == spnames[x], ], ]$ID
      sps[[x]] <- na.omit(data.frame(ID, Species = spnames[x]))
    }
    close(pb)

    sps <- do.call(rbind, sps)
  }

  return(sps)
}



#' Creates a data.frame of species' references from a list of raster layers
#'
#' @description Creates a data.frame of species' references that contains
#' longitude, latitude, and species name, using a list of raster layers as
#' input. Useful when raster layers have distinct extent or resolution.
#'
#' @param raster_list list of SpatRaster objects. Each raster layer must be
#' named as the species that it represents, and values in each layer must be
#' 1 (presence) and 0 (absence).
#'
#' @return
#' A data.frame of species geographic records derived from values of presence
#' in each layer from the list of raster layers.
#'
#' @usage
#' rlist_2data(raster_list)
#'
#' @export
#' @importFrom terra as.data.frame
#'
#' @examples
#' # Data
#' rsp <- terra::rast(system.file("extdata/sp_layers.tif",
#'                                package = "biosurvey"))
#' names(rsp) <- paste0("Species_", 1:5)
#'
#' rlist <- lapply(1:5, function(x) {rsp[[x]]})
#'
#' # Species data from RasterStack
#' sp_data <- rlist_2data(raster_list = rlist)
#' summary(sp_data)

rlist_2data <- function(raster_list) {
  # Initial tests
  if (missing(raster_list)) {
    stop("Argument 'raster_list' must be defined")
  }
  if (!is.list(raster_list)) {
    stop("'raster_list' must be a list of raster layers")
  }
  inclas <- sapply(raster_list, function(x) {class(x)[1] != "SpatRaster"})
  if (any(inclas)) {
    stop("All elements in 'raster_list' must be of class 'SpatRaster'")
  }

  # Running in loop for all elements of list
  sps <- lapply(raster_list, function(x) {
    # Raster to data.frame
    sppm <- terra::as.data.frame(x, xy = TRUE)

    # Preparing data
    data.frame(sppm[sppm[, 3] == 1, 1:2], names(x))
  })

  sps <- do.call(rbind, sps)

  colnames(sps) <- c("Longitude", "Latitude", "Species")

  return(sps)
}



#' Creates a data.frame of species' references from files in a directory
#'
#' @description Creates a data.frame of species' references that contains
#' longitude, latitude, and species name, from a character.
#'
#' @param path (character) full path name of directory containing raster,
#' shapefiles or geopackage files representing species geographic
#' ranges. Each file must be named as the species that it represents. All files
#' must be in the same format. If files are raster layers, values in each layer
#' must be 1 (presence, suitable) and 0 (absence, unsuitable).
#' @param format (character) the format files found in \code{path}. Current
#' available formats are: "shp", "gpkg", "GTiff", and "ascii".
#' @param spdf_grid geographic grid for the region of interest (output of
#' function \code{\link{grid_from_region}}). Used when format equals "shp" or
#' "gpkg". Default = NULL.
#' @param parallel (logical) whether to perform analyses in parallel.
#' Default = FALSE.
#' @param n_cores (numeric) number of cores to be used when \code{parallel} =
#' TRUE. The default, NULL, uses available cores - 1.
#'
#' @return
#' If files are in raster format, a data.frame of species geographic records
#' derived from values of presence in each layer.
#'
#' If files are not in raster format, a data.frame of species' found in distinct
#' positions (defined with identifiers); includes two columns: "ID" and
#' "Species".
#'
#' @usage
#' files_2data(path, format, spdf_grid = NULL, parallel = FALSE, n_cores = NULL)
#'
#' @export
#' @importFrom terra vect
#' @importFrom sp over
#' @importFrom terra rast as.data.frame
#' @importFrom stats na.omit
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @examples
#' \donttest{
#' # Data for examples
#' data("mx", package = "biosurvey")
#' data("species_data", package = "biosurvey")
#'
#' # Saving species data in a temporal directory
#' tdir <- file.path(tempdir(), "testbio")
#' dir.create(tdir)
#'
#' namessp <- paste0("species_", 1:length(species_data))
#'
#'
#' for (i in 1:length(species_data)) {
#'   rgdal::writeOGR(species_data[i, ], dsn = tdir, layer = namessp[i],
#'                   driver = "ESRI Shapefile")
#' }
#'
#' # Preparing grid for analysis
#' grid_reg <- grid_from_region(region = mx, cell_size = 100)
#'
#' # Running analysis with data from directory
#' sp_data <- files_2data(path = tdir, format = "shp", spdf_grid = grid_reg)
#' }

files_2data <- function(path, format, spdf_grid = NULL, parallel = FALSE,
                        n_cores = NULL) {
  # Initial tests
  if (missing(path)) {
    stop("Argument 'path' must be defined")
  }
  if (missing(format)) {
    stop("Argument 'format' must be defined")
  }
  if (!format %in% c("shp", "gpkg", "GTiff", "ascii")) {
    stop(paste("'format'", format, "is not supported, see function's help"))
  }

  if (is.null(spdf_grid)) {
    stop("Argument 'spdf_grid' must be defined if 'format' is 'shp' or 'gpkg'")
  }

  # spdf_grid crs
  crsg <- terra::crs(spdf_grid)

  # Finding files according to format
  if (format %in% c("GTiff", "ascii")) {
    format <- match_rformat(format)
  }
  patt <- paste0(".", format, "$")

  mlist <- list.files(path = path, pattern = patt, full.names = TRUE)
  spnames <- gsub(patt, "", list.files(path = path, pattern = patt))


  if (length(mlist) == 0) {
    stop(paste("No files were found in", path, "with the extension specified in 'format'"))
  }

  if (parallel == TRUE) {
    ## Preparing parallel running
    n_cores <- ifelse(is.null(n_cores), parallel::detectCores() - 1, n_cores)

    ## Progress combine (rbind) function
    pb <- utils::txtProgressBar(min = 1, max = length(spnames), style = 3)
    progress <- function(n) {
      utils::setTxtProgressBar(pb, n)
    }
    opts <- list(progress = progress)

    ## Make cluster
    cl <- snow::makeSOCKcluster(n_cores)
    doSNOW::registerDoSNOW(cl)

    ## wrap vectors
    spdf_grid <- terra::wrap(spdf_grid)

    ## Processing
    sps <- foreach::foreach(
      i = 1:length(spnames), .inorder = FALSE,
      .combine = "rbind", .options.snow = opts
    ) %dopar% {
      if (format %in% c("shp", "gpkg")) {
        rs <- terra::vect(mlist[i])

        if (terra::crs(rs) != crsg) {
          rs <- terra::project(rs, crsg)
        }

        ID <- terra::unwrap(spdf_grid)[rs, ]$ID

        if (length(ID) > 0) {
          return(na.omit(data.frame(ID = ID, Species = spnames[i])))
        } else {
          return(na.omit(data.frame(ID = NA, Species = NA)))
        }

      } else {
        rs <- terra::rast(mlist[x])

        if (terra::crs(rs) != crsg) {
          rs <- terra::project(rs, crsg)
        }

        # Raster to data.frame
        sppm <- terra::as.data.frame(rs, xy = TRUE)

        # Preparing data
        return(data.frame(sppm[sppm[, 3] == 1, 1:2], spnames[x]))
      }
    }
    snow::stopCluster(cl)

  } else {
    ## Progress bar
    pb <- utils::txtProgressBar(min = 1, max = length(spnames), style = 3)

    # Running in loop for all elements of list
    sps <- list()

    for (x in 1:length(spnames)) {
      Sys.sleep(0.1)
      utils::setTxtProgressBar(pb, x)

      if (format %in% c("shp", "gpkg")) {
        rs <- terra::vect(mlist[x])

        if (terra::crs(rs) != crsg) {
          rs <- terra::project(rs, crsg)
        }

        ID <- spdf_grid[rs, ]$ID

        if (length(ID) > 0) {
          sps[[x]] <- na.omit(data.frame(ID, Species = spnames[x]))
        } else {
          sps[[x]] <- na.omit(data.frame(ID = NA, Species = NA))
        }

      } else {
        rs <- terra::rast(mlist[x])

        if (terra::crs(rs) != crsg) {
          rs <- terra::project(rs, crsg)
        }

        # Raster to data.frame
        sppm <- terra::as.data.frame(rs, xy = TRUE)

        # Preparing data
        sps[[x]] <- data.frame(sppm[sppm[, 3] == 1, 1:2], spnames[x])
      }
    }
    close(pb)

    sps <- do.call(rbind, sps)
  }

  if (format %in% c("shp", "gpkg")) {
    colnames(sps) <- c("ID", "Species")
  } else {
    colnames(sps) <- c("Longitude", "Latitude", "Species")
  }

  return(sps)
}



#' Creates presence-absence matrix from a data.frame
#'
#' @description Creates a presence-absence matrix (PAM) from a data.frame that
#' contains species names and identifiers of positions where species are found.
#'
#' @param data data.frame of species found in distinct positions (defined by
#' identifiers). Must include at least two columns: "ID" and "Species".
#' @param ID_column (character) name of the column containing identifiers.
#' @param species_column (character) name of the column containing species
#' names.
#'
#' @return
#' Species' presence (1) and absence (0) matrix for a set of positions defined
#' by identifiers.
#'
#' @usage
#' PAM_from_table(data, ID_column, species_column)
#'
#' @export
#'
#' @examples
#' # Data
#' data("sp_data", package = "biosurvey")
#'
#' # PAM
#' pam <- PAM_from_table(data = sp_data, ID_column = "ID",
#'                       species_column = "Species")
#' pam[1:10, c(1, 21:25)]

PAM_from_table <- function(data, ID_column, species_column) {
  # Initial tests
  if (missing(data)) {
    stop("Argument 'data' must be defined")
  }
  if (missing(ID_column)) {
    stop("Argument 'ID_column' must be defined")
  }
  if (missing(species_column)) {
    stop("Argument 'species_column' must be defined")
  }
  if(class(data[, ID_column])[1] != "factor") {
    data[, ID_column] <- as.factor(data[, ID_column])
  }

  # Transform species column to characters
  allsp <- as.character(unique(data[, species_column]))

  # Count of the species per ID
  counts <- sapply(allsp, function(x) {
    spp <- data[as.character(data[, species_column]) == x, ]
    table(spp[, ID_column])
  })

  # Fixed details of PAM
  counts[counts > 0] <- 1
  nams <- colnames(counts)
  counts <- data.frame(rownames(counts), counts)
  colnames(counts) <- c(ID_column, nams)

  return(counts)
}





#' Helper to subset PAM according to selected sites
#' @param selected_sites list of selected sites. See any of the main functions
#' to perform such a selection: \code{\link{random_selection}},
#' \code{\link{uniformG_selection}}, \code{\link{uniformE_selection}}, or
#' \code{\link{EG_selection}}.
#' @param base_PAM object of class base_PAM obtained using the function
#' \code{\link{prepare_base_PAM}}.
#'
#' @return
#' A list of selected site data.frames with information of PAM added as
#' additional columns.
#'
#' @export
#' @importFrom terra crs vect as.data.frame

selected_sites_PAM <- function(selected_sites, base_PAM) {
  # Initial tests
  if(missing(selected_sites)) {
    stop("Argument 'selected_sites' must be defined.")
  }
  if (missing(base_PAM)) {
    stop("Argument 'base_PAM' must be defined.")
  }

  WGS84 <- terra::crs(base_PAM$PAM)

  # Matching sites with PAM IDs
  ls <- lapply(selected_sites, function(x) {
    xp <- terra::vect(x, geom = c("Longitude", "Latitude"), crs = WGS84)
    xid <- data.frame(ID = base_PAM$PAM[xp, ]$ID, x)
    pam <- terra::as.data.frame(terra::vect(base_PAM$PAM))
    colnames(pam)[2:3] <- c("Longitude_PAM", "Latitude_PAM")
    colnames(xid)[2:3] <- c("Longitude_master", "Latitude_master")
    merge(xid, pam, by = "ID")
  })
  names(ls) <- names(selected_sites)
  return(ls)
}




#' Helper to refill a list of PAM indices with new or more results
#'
#' @param initial_index_list list of PAM indices to be refill. Indices present
#' in this list and absent in \code{new_index_list} are maintained.
#' @param new_index_list list of PAM indices to be used to refill
#' \code{initial_index_list}. New indices are included in the resulting list.
#' Indices present in both lists are updated using the values of this list.
#'
#' @export
#'
#' @return
#' A list of PAM indices containing old and new values for its indices.

refill_PAM_indices <- function(initial_index_list, new_index_list) {
  # Initial test
  if (missing(initial_index_list)) {
    stop("Argument 'initial_index_list' must be defined.")
  }
  if (missing(new_index_list)) {
    stop("Argument 'new_index_list' must be defined.")
  }

  # Starting filling list
  index_list <- initial_index_list

  # Non basic
  ## One value
  index_list$One_value_indices["Av_dispersion_field", ] <- ifelse(
    is.na(new_index_list$One_value_indices["Av_dispersion_field", ]) &
      !is.na(initial_index_list$One_value_indices["Av_dispersion_field", ]),
    initial_index_list$One_value_indices["Av_dispersion_field", ],
    new_index_list$One_value_indices["Av_dispersion_field", ]
  )

  index_list$One_value_indices["Av_diversity_field", ] <- ifelse(
    is.na(new_index_list$One_value_indices["Av_diversity_field", ]) &
      !is.na(initial_index_list$One_value_indices["Av_diversity_field", ]),
    initial_index_list$One_value_indices["Av_diversity_field", ],
    new_index_list$One_value_indices["Av_diversity_field", ]
  )

  index_list$One_value_indices["Av_shared_community_composition", ] <- ifelse(
    is.na(new_index_list$One_value_indices["Av_shared_community_composition", ]) &
      !is.na(initial_index_list$One_value_indices["Av_shared_community_composition", ]),
    initial_index_list$One_value_indices["Av_shared_community_composition", ],
    new_index_list$One_value_indices["Av_shared_community_composition", ]
  )

  index_list$One_value_indices["Additive_Beta", ] <- ifelse(
    is.na(new_index_list$One_value_indices["Additive_Beta", ]) &
      !is.na(initial_index_list$One_value_indices["Additive_Beta", ]),
    initial_index_list$One_value_indices["Additive_Beta", ],
    new_index_list$One_value_indices["Additive_Beta", ]
  )

  index_list$One_value_indices["Beta_Whittaker", ] <- ifelse(
    is.na(new_index_list$One_value_indices["Beta_Whittaker", ]) &
      !is.na(initial_index_list$One_value_indices["Beta_Whittaker", ]),
    initial_index_list$One_value_indices["Beta_Whittaker", ],
    new_index_list$One_value_indices["Beta_Whittaker", ]
  )

  index_list$One_value_indices["Beta_Legendre", ] <- ifelse(
    is.na(new_index_list$One_value_indices["Beta_Legendre", ]) &
      !is.na(initial_index_list$One_value_indices["Beta_Legendre", ]),
    initial_index_list$One_value_indices["Beta_Legendre", ],
    new_index_list$One_value_indices["Beta_Legendre", ]
  )

  index_list$One_value_indices["Schluter_cov_sites_composition", ] <- ifelse(
    is.na(new_index_list$One_value_indices["Schluter_cov_sites_composition", ]) &
      !is.na(initial_index_list$One_value_indices["Schluter_cov_sites_composition", ]),
    initial_index_list$One_value_indices["Schluter_cov_sites_composition", ],
    new_index_list$One_value_indices["Schluter_cov_sites_composition", ]
  )

  index_list$One_value_indices["Schluter_cov_species_ranges", ] <- ifelse(
    is.na(new_index_list$One_value_indices["Schluter_cov_species_ranges", ]) &
      !is.na(initial_index_list$One_value_indices["Schluter_cov_species_ranges", ]),
    initial_index_list$One_value_indices["Schluter_cov_species_ranges", ],
    new_index_list$One_value_indices["Schluter_cov_species_ranges", ]
  )
  index_list$One_value_indices["Wright_Reeves_nestedness", ] <- ifelse(
    is.na(new_index_list$One_value_indices["Wright_Reeves_nestedness", ]) &
      !is.na(initial_index_list$One_value_indices["Wright_Reeves_nestedness", ]),
    initial_index_list$One_value_indices["Wright_Reeves_nestedness", ],
    new_index_list$One_value_indices["Wright_Reeves_nestedness", ]
  )
  index_list$One_value_indices["Stone_Roberts_Cscore", ] <- ifelse(
    is.na(new_index_list$One_value_indices["Stone_Roberts_Cscore", ]) &
      !is.na(initial_index_list$One_value_indices["Stone_Roberts_Cscore", ]),
    initial_index_list$One_value_indices["Stone_Roberts_Cscore", ],
    new_index_list$One_value_indices["Stone_Roberts_Cscore", ]
  )

  ## Lists
  if (all(is.na(new_index_list$Dispersion_field)) &
      any(!is.na(initial_index_list$Dispersion_field))) {
    index_list$Dispersion_field <- initial_index_list$Dispersion_field
  } else {
    index_list$Dispersion_field <- new_index_list$Dispersion_field
  }

  if (all(is.na(new_index_list$Diversity_field)) &
      any(!is.na(initial_index_list$Diversity_field))) {
    index_list$Diversity_field <- initial_index_list$Diversity_field
  } else {
    index_list$Diversity_field <- new_index_list$Diversity_field
  }

  if (all(is.na(new_index_list$Shared_community_composition)) &
      any(!is.na(initial_index_list$Shared_community_composition))) {
    index_list$Shared_community_composition <- initial_index_list$Shared_community_composition
  } else {
    index_list$Shared_community_composition <- new_index_list$Shared_community_composition
  }

  if (all(is.na(new_index_list$Mean_composition_covariance)) &
      any(!is.na(initial_index_list$Mean_composition_covariance))) {
    index_list$Mean_composition_covariance <- initial_index_list$Mean_composition_covariance
  } else {
    index_list$Mean_composition_covariance <- new_index_list$Mean_composition_covariance
  }

  if (all(is.na(new_index_list$Mean_range_covariance)) &
      any(!is.na(initial_index_list$Mean_range_covariance))) {
    index_list$Mean_range_covariance <- initial_index_list$Mean_range_covariance
  } else {
    index_list$Mean_range_covariance <- new_index_list$Mean_range_covariance
  }

  if (all(is.na(new_index_list$Cov_mat_sites_composition)) &
      any(!is.na(initial_index_list$Cov_mat_sites_composition))) {
    index_list$Cov_mat_sites_composition <- initial_index_list$Cov_mat_sites_composition
  } else {
    index_list$Cov_mat_sites_composition <- new_index_list$Cov_mat_sites_composition
  }

  if (all(is.na(new_index_list$Cov_mat_species_ranges)) &
      any(!is.na(initial_index_list$Cov_mat_species_ranges))) {
    index_list$Cov_mat_species_ranges <- initial_index_list$Cov_mat_species_ranges
  } else {
    index_list$Cov_mat_species_ranges <- new_index_list$Cov_mat_species_ranges
  }

  return(index_list)
}
