#' Creates grid for a given geographic region
#'
#' @description divides the region of interest in a grid of a specific cell size.
#'
#' @param region SpatialPolygonsDataFrame of the region of interest. Object must
#' be unprojected, World Geodetic System (WGS84).
#' @param cell_size (numeric) resolution for grid (single number or vector of
#' two numbers) in kilometers (km).
#' @param complete_cover (logical) whether or not to include cells of grid
#' partially overlapped with region. Default = TRUE.
#'
#' @return
#' Gridded SpatialPolygonsDataFrame for the region of interest. Each grid cell
#' is related to a specific ID and longitude and latitude coordinates.
#'
#' @usage
#' grid_from_region(region, cell_size, complete_cover = TRUE)
#'
#' @export
#' @importFrom raster extent raster res values mask rasterToPolygons rasterToPoints
#' @importFrom raster projectRaster rasterize
#' @importFrom sp proj4string CRS spTransform
#' @importFrom rgeos gCentroid
#'
#' @examples
#' # Data
#' data("mx", package = "biosurvey")
#'
#' # Create grid from polygon
#' grid_reg <- grid_from_region(region = mx, cell_size = 100)
#'
#' sp::plot(grid_reg)
#' grid_reg

grid_from_region <- function(region, cell_size, complete_cover = TRUE) {
  # Initial tests
  if (missing(region)) {
    stop("Argument 'region' must be defined")
  }
  if (missing(cell_size)) {
    stop("Argument 'cell_size' must be defined")
  } else {

    # Projecting region toLambert equeal area projection
    if (is.na(sp::proj4string(region))) {
      stop("'region' must be projected to WGS84 (EPSG:4326)")
    }
    WGS84 <- sp::CRS("+init=epsg:4326")
    region <- sp::spTransform(region, WGS84)
    cent <- rgeos::gCentroid(region, byid = FALSE)@coords
    LAEA <- sp::CRS(paste0("+proj=laea +lat_0=", cent[2], " +lon_0=", cent[1],
                           " +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
    region <- sp::spTransform(region, LAEA)

    # Test if dimensions are valid
    dims <- raster::extent(region)
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
  grid <- raster::raster(raster::extent(region))

  # Grid resolution and values
  raster::res(grid) <- cell_size * 1000
  raster::values(grid) <- 1

  # Grid projection
  sp::proj4string(grid) <- sp::proj4string(region)

  # Extract grid with region
  if (complete_cover == TRUE) {
    SpP_ras <- raster::rasterize(region, grid, getCover = TRUE)
    SpP_ras[SpP_ras == 0] <- NA
    grid_reg <- raster::mask(grid, SpP_ras)
  } else {
    message("Cells partially covered by polygon representing region won't be included.",
            "\nTo include such cells use 'complete_cover' = TRUE.")
    grid_reg <- raster::mask(grid, region)
  }

  # Back to WGS84
  grid_reg <- raster::projectRaster(grid_reg, crs = WGS84)

  # Grid for region of interest
  grid_r_pol <- raster::rasterToPolygons(grid_reg)

  # Points for region of interest
  matrix_a <- raster::rasterToPoints(grid_reg)

  # Adding ID for PAM
  ID <- raster::extract(grid_reg, matrix_a[, 1:2], cellnumbers = TRUE)[, 1]
  grid_r_pol@data <- data.frame(ID = ID, Longitude = matrix_a[, 1],
                                Latitude = matrix_a[, 2])

  return(grid_r_pol)
}



#' Creates a data.frame of species' references from RasterStack
#'
#' @description creates a data.frame of species' references that contains
#' longitude, latitude, and species name, using a RasterStack or a RasterBrick
#' as input.
#'
#' @param species_layers RasterStack or RasterBrick object. Each layer must be
#' named as the species that it represents, and values in each layer must be
#' 1 (presence) and 0 (absence).
#'
#' @return
#' A data.frame of species geographic records derived from values of presence
#' in each layer from the RasterStack.
#'
#' @usage
#' stack_2data(species_layers)
#'
#' @export
#' @importFrom raster rasterToPoints
#'
#' @examples
#' # Data
#' rsp <- raster::stack(system.file("extdata/sp_layers.tif",
#'                                  package = "biosurvey"))
#' names(rsp) <- paste0("Species_", 1:5)
#'
#' # Species data from RasterStack
#' sp_data <- stack_2data(species_layers = rsp)
#' summary(sp_data)

stack_2data <- function(species_layers) {
  # Initial tests
  if (missing(species_layers)) {
    stop("Argument 'species_layers' must be defined")
  }
  if (class(species_layers)[1] != "RasterStack") {
    stop("'species_layers' must be of class 'RasterStack'")
  }

  # Stack to matrix
  sppm <- raster::rasterToPoints(species_layers)
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



#' Creates a data.frame of species' references from SpatialPolygonsDataFrame
#'
#' @description creates a data.frame of species' references that contains
#' identifiers of portion and species name, using a SpatialPolygonsDataFrame as
#' input.
#'
#' @param spdf_object SpatialPolygonsDataFrame representing species' geographic
#' distributions. The data.frame associated with the object must contain a
#' column named "Species" to distinguish among features.
#' @param spdf_grid geographic grid for the region of interest (output of
#' function \code{\link{grid_from_region}}).
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
#' @importFrom utils txtProgressBar setTxtProgressBar flush.console
#' @importFrom sp over
#' @importFrom stats na.omit
#'
#' @examples
#' # Data
#' data("species_data", package = "biosurvey")
#' data("mx", package = "biosurvey")
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
  cond <- c(class(spdf_object)[1] != "SpatialPolygonsDataFrame",
            class(spdf_grid)[1] != "SpatialPolygonsDataFrame")
  if (any(cond)) {
    stop("'spdf_object' and 'spdf_grid' must be of class 'SpatialPolygonsDataFrame'")
  }

  # Names to be matched
  ID <- spdf_grid@data$ID
  spnames <- as.character(spdf_object@data$Species)

  if (parallel == TRUE) {
    ## Preparing parallel running
    n_cores <- ifelse(is.null(n_cores), parallel::detectCores() - 1, n_cores)

    ## Progress combine (rbind) function
    fpc <- function(iterator){
      pb <- utils::txtProgressBar(min = 1, max = iterator - 1, style = 3)
      count <- 0
      function(...) {
        count <<- count + length(list(...)) - 1
        utils::setTxtProgressBar(pb, count)
        utils::flush.console()
        Sys.sleep(0.1)
        rbind(...)
      }
    }

    ## Start a cluster
    cl <- parallel::makeCluster(n_cores, type = 'SOCK')
    doParallel::registerDoParallel(cl)

    ## Processing
    sps <- foreach::foreach(i = 1:length(spnames), .inorder = TRUE,
                            .combine = fpc(length(spnames))) %dopar% {
                              sp <- sp::over(spdf_grid,
                                             spdf_object[spnames == spnames[i], ])
                              return(na.omit(data.frame(ID, sp)))
                            }

    parallel::stopCluster(cl)
  } else {
    ## Progress bar
    pb <- utils::txtProgressBar(min = 1, max = length(spnames), style = 3)

    # Running in loop for all elements of list
    sps <- list()

    for (x in 1:length(spnames)) {
      Sys.sleep(0.1)
      utils::setTxtProgressBar(pb, x)

      # Preparing data
      sp <- sp::over(spdf_grid, spdf_object[spnames == spnames[x], ])
      sps[[x]] <- na.omit(data.frame(ID, sp))
    }

    sps <- do.call(rbind, sps)
  }

  colnames(sps) <- c("ID", "Species")

  return(sps)
}



#' Creates a data.frame of species' references from a list of raster layers
#'
#' @description creates a data.frame of species' references that contains
#' longitude, latitude, and species name, using a list of raster layers as
#' input. Useful when raster layers have distinct extent or resolution.
#'
#' @param raster_list list of RasterLayer objects. Each raster layer must be
#' named as the species that it represents, and values in each layer must be
#' 1 (presence) and 0 (absence).
#' @param parallel (logical) whether to perform analyses in parallel.
#' Default = FALSE.
#' @param n_cores (numeric) number of cores to be used when \code{parallel} =
#' TRUE. The default, NULL, uses available cores - 1.
#'
#' @return
#' A data.frame of species geographic records derived from values of presence
#' in each layer from the list of raster layers.
#'
#' @usage
#' rlist_2data(raster_list, parallel = FALSE, n_cores = NULL)
#'
#' @export
#' @importFrom raster rasterToPoints
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom utils txtProgressBar setTxtProgressBar flush.console
#'
#' @examples
#' # Data
#' rsp <- raster::stack(system.file("extdata/sp_layers.tif",
#'                      package = "biosurvey"))
#' names(rsp) <- paste0("Species_", 1:5)
#'
#' rlist <- lapply(1:5, function(x) {rsp[[x]]})
#'
#' # Species data from RasterStack
#' sp_data <- rlist_2data(raster_list = rlist)
#' summary(sp_data)

rlist_2data <- function(raster_list, parallel = FALSE, n_cores = NULL) {
  # Initial tests
  if (missing(raster_list)) {
    stop("Argument 'raster_list' must be defined")
  }
  if (!is.list(raster_list)) {
    stop("'raster_list' must be a list of raster layers")
  }
  inclas <- sapply(raster_list, function(x) {class(x)[1] != "RasterLayer"})
  if (any(inclas)) {
    stop("All elements in 'raster_list' must be of class 'RasterLayer'")
  }

  # Running in loop for all elements of list
  sps <- lapply(1:length(raster_list), function(x) {
    # Raster to matrix
    sppm <- raster::rasterToPoints(raster_list[[x]])
    spname <- names(raster_list[[x]])

    # Preparing data
    data.frame(sppm[sppm[, 3] == 1, 1:2], spname)
  })

  if (parallel == TRUE) {
    ## Preparing parallel running
    n_cores <- ifelse(is.null(n_cores), parallel::detectCores() - 1, n_cores)

    ## Progress combine (rbind) function
    fpc <- function(iterator){
      pb <- utils::txtProgressBar(min = 1, max = iterator - 1, style = 3)
      count <- 0
      function(...) {
        count <<- count + length(list(...)) - 1
        utils::setTxtProgressBar(pb, count)
        utils::flush.console()
        Sys.sleep(0.1)
        rbind(...)
      }
    }

    ## Start a cluster
    cl <- parallel::makeCluster(n_cores, type = 'SOCK')
    doParallel::registerDoParallel(cl)

    ## Processing
    sps <- foreach::foreach(i = 1:length(raster_list), .inorder = TRUE,
                            .combine = fpc(length(raster_list))) %dopar% {
                              # Raster to matrix
                              sppm <- raster::rasterToPoints(raster_list[[i]])
                              spname <- names(raster_list[[i]])

                              # Preparing data
                              cond <- sppm[, 3] == 1
                              data.frame(sppm[cond, 1], sppm[cond, 2], spname)
                            }

    parallel::stopCluster(cl)
  } else {
    ## Progress bar
    pb <- utils::txtProgressBar(min = 1, max = length(raster_list), style = 3)

    # Running in loop for all elements of list
    sps <- list()

    for (x in 1:length(raster_list)) {
      Sys.sleep(0.1)
      utils::setTxtProgressBar(pb, x)

      # Raster to matrix
      sppm <- raster::rasterToPoints(raster_list[[x]])
      spname <- names(raster_list[[x]])

      # Preparing data
      cond <- sppm[, 3] == 1
      sps[[x]] <- data.frame(sppm[cond, 1], sppm[cond, 2], spname)
    }

    sps <- do.call(rbind, sps)
  }

  colnames(sps) <- c("Longitude", "Latitude", "Species")

  return(sps)
}



#' Creates a data.frame of species' references from files in a directory
#'
#' @description creates a data.frame of species' references that contains
#' longitude, latitude, and species name, from a character.
#'
#' @param path (character) full path name of directory containing raster,
#' shapefiles, geopackage, or GeoJSON files representing species geographic
#' ranges. Each file must be named as the species that it represents. All files
#' must be in the same format. If files are raster, values in each layer must be
#' 1 (presence) and 0 (absence).
#' @param format (character) the format files found in \code{path}. Current
#' available formats are: "shp", "gpkg", "geojson", "GTiff", and "ascii".
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
#' @importFrom rgdal readOGR
#' @importFrom sp over
#' @importFrom raster raster rasterToPoints
#' @importFrom stats na.omit
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom utils txtProgressBar setTxtProgressBar flush.console
#'
#' @examples
#' # example of how to define arguments, check argument descriptions above
#'
#' \dontrun{
#' # Using folder with rasters in GeoTiff format
#' sp_data <- files_2data(path = "Folder_with_rasters", format = "GTiff")
#' summary(sp_data)
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
  if (!format %in% c("shp", "gpkg", "geojson", "GTiff", "ascii")) {
    stop(paste("'format'", format, "is not supported, see function's help"))
  }

  # Finding files according to format
  if (format %in% c("shp", "gpkg", "geojson")) {
    if (is.null(spdf_grid)) {
      stop("Argument 'spdf_grid' must be defined if 'format' is shp, gpkg, or geojson")
    }
    # Names to be matched
    ID <- spdf_grid@data$ID

    if (format == "shp") {
      patt <- ".shp$"
      subs <- ".shp"
      mlist <- gsub(subs, "", list.files(path = path, pattern = patt))
      spnames <- mlist
    } else {
      if (format == "gpkg") {
        patt <- ".gpkg$"
        subs <- ".gpkg"
        mlist <- list.files(path = path, pattern = patt)
        spnames <- gsub(subs, "", mlist)
      } else {
        patt <- ".geojson$"
        subs <- ".geojson"
        mlist <- list.files(path = path, pattern = patt)
        spnames <- gsub(subs, "", mlist)
      }
    }
  } else {
    subs <- match_rformat(format)
    patt <- paste0(subs, "$")
    mlist <- list.files(path = path, pattern = patt, full.names = TRUE)
    spnames <- gsub(subs, "", list.files(path = path, pattern = patt))
  }

  if (length(mlist) == 0) {
    stop(paste("No file was found in", path, "with the extension specified in 'format'"))
  }

  if (parallel == TRUE) {
    ## Preparing parallel running
    n_cores <- ifelse(is.null(n_cores), parallel::detectCores() - 1, n_cores)

    ## Progress combine (rbind) function
    fpc <- function(iterator){
      pb <- utils::txtProgressBar(min = 1, max = iterator - 1, style = 3)
      count <- 0
      function(...) {
        count <<- count + length(list(...)) - 1
        utils::setTxtProgressBar(pb, count)
        utils::flush.console()
        Sys.sleep(0.1)
        rbind(...)
      }
    }

    ## Start a cluster
    cl <- parallel::makeCluster(n_cores, type = 'SOCK')
    doParallel::registerDoParallel(cl)

    ## Processing
    sps <- foreach::foreach(i = 1:length(spnames), .inorder = TRUE,
                            .combine = fpc(length(spnames))) %dopar% {
                              if (format %in% c("shp", "gpkg", "geojson")) {
                                ## Reading data
                                if (format == "shp") {
                                  rs <- rgdal::readOGR(dsn = path,
                                                       layer = mlist[i],
                                                       verbose = FALSE)
                                } else {
                                  if (format == "gpkg") {
                                    rs <- rgdal::readOGR(paste0(path, "/",
                                                                mlist[i]),
                                                         spnames[x],
                                                         verbose = FALSE)
                                  } else {
                                    rs <- rgdal::readOGR(paste0(path, "/",
                                                                mlist[i]),
                                                         verbose = FALSE)
                                  }
                                }

                                ## Preparing data
                                sp <- sp::over(spdf_grid, rs)[, 1]
                                sppm <- na.omit(data.frame(ID, Species = sp))
                                sppm$Species <- spnames[i]
                                return(sppm)

                              } else {
                                ## Raster from file
                                rs <- raster::raster(mlist[i])

                                ## Raster to matrix
                                sppm <- raster::rasterToPoints(rs)

                                ## Preparing data
                                cond <- sppm[, 3] == 1
                                return(data.frame(sppm[cond, 1], sppm[cond, 2],
                                                  spnames[i]))
                              }
                            }

    parallel::stopCluster(cl)
  } else {
    ## Progress bar
    pb <- utils::txtProgressBar(min = 1, max = length(spnames), style = 3)

    # Running in loop for all elements of list
    sps <- list()

    for (x in 1:length(spnames)) {
      Sys.sleep(0.1)
      utils::setTxtProgressBar(pb, x)

      if (format %in% c("shp", "gpkg", "geojson")) {
        ## Reading data
        if (format == "shp") {
          rs <- rgdal::readOGR(dsn = path, layer = mlist[x], verbose = FALSE)
        } else {
          if (format == "gpkg") {
            rs <- rgdal::readOGR(paste0(path, "/", mlist[x]), spnames[x],
                                 verbose = FALSE)
          } else {
            rs <- rgdal::readOGR(paste0(path, "/", mlist[x]),
                                 verbose = FALSE)
          }
        }

        ## Preparing data
        sppm <- na.omit(data.frame(ID, Species = sp::over(spdf_grid, rs)[, 1]))
        sppm$Species <- spnames[x]
        sps[[x]] <- sppm

      } else {
        ## Raster from file
        rs <- raster::raster(mlist[x])

        ## Raster to matrix
        sppm <- raster::rasterToPoints(rs)

        ## Preparing data
        cond <- sppm[, 3] == 1
        sps[[x]] <- data.frame(sppm[cond, 1], sppm[cond, 2], spnames[x])
      }
    }

    sps <- do.call(rbind, sps)
  }

  if (format %in% c("shp", "gpkg", "geojson")) {
    colnames(sps) <- c("ID", "Species")
  } else {
    colnames(sps) <- c("Longitude", "Latitude", "Species")
  }

  return(sps)
}



#' Creates presence-absence matrix from a data.frame
#'
#' @description creates a presence-absence matrix (PAM) from a data.frame that
#' contains species names and identifiers of positions where species are found.
#'
#' @param data data.frame of species' found in distinct positions (defined by
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
#' @importFrom sp CRS over SpatialPointsDataFrame

selected_sites_PAM <- function(selected_sites, base_PAM) {
  # Initial tests
  if(missing(selected_sites)) {
    stop("Argument 'selected_sites' must be defined.")
  }
  if (missing(base_PAM)) {
    stop("Argument 'base_PAM' must be defined.")
  }

  WGS84 <- sp::CRS("+init=epsg:4326")

  # Matching sites with PAM IDs
  ls <- lapply(selected_sites, function(x) {
    xp <- sp::SpatialPointsDataFrame(x[, 1:2], x, proj4string = WGS84)
    xid <- data.frame(sp::over(xp, base_PAM$PAM[, "ID"]), x)
    pam <- base_PAM$PAM@data
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
