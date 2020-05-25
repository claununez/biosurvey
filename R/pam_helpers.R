#' Creates the grids in the geography to be used to extract the values of the base_pam
#'
#' @description Divides the region of interest in grids of a determined size that
#' will be used as the cells to extract the values of the base_pam function.
#'
#' @param region SpatialPolygonsDataFrame of the region of interest.
#' @param cell_size Resolution for grid (single number or vector of two numbers)
#' in decimal degrees.
#'
#' @return
#'
#' @usage
#' grid_from_region(region, cell_size)
#'
#' @export
#' @importFrom raster extent raster res values mask rasterToPolygons rasterToPoints
#' @importFrom sp proj4string

grid_from_region <- function(region, cell_size) {
  # Initial tests
  if (missing(region)) {
    stop("Argument 'region' must be defined")
  }
  if (missing(cell_size)) {
    stop("Argument 'cell_size' must be defined")
  }
  dims <- raster::extent(region)
  xdim <- diff(dims[1:2])
  ydim <- diff(dims[3:4])
  if (cell_size >= xdim & cell_size >= ydim) {
    stop("'cell_size' must be smaller than at least one of the dimensions of 'region'")
  }
  if (is.na(sp::proj4string(region))) {
    stop("'region' must be projected to WGS84 (EPSG:4326)")
  }

  # creating a grid
  grid <- raster::raster(raster::extent(region))

  # grid resolution and values
  raster::res(grid) <- cell_size
  raster::values(grid) <- 0

  # grid projection
  sp::proj4string(grid) <- sp::proj4string(region)

  # extract grid with region
  grid_reg <- raster::mask(grid, region)

  # grid for region of interest
  grid_r_pol <- raster::rasterToPolygons(grid_reg)

  # points for region of interest
  matrix_a <- raster::rasterToPoints(grid_reg)

  # Adding ID for PAM
  ID <- raster::extract(grid_reg, matrix_a[, 1:2], cellnumbers = TRUE)[, 1]
  grid_r_pol@data <- data.frame(ID = ID, Longitude = matrix_a[, 1],
                                Latitude = matrix_a[, 2])
  return(grid_r_pol)
}





#' Creates PAM from a data frame of species references
#'
#' @description Creates a presence-absence matrix (PAM) from a data frame that
#' contains identifiers and species names.
#'
#' @param data data frame objects of species' presence (1) absence (0)
#' @param ID_column
#' @param species_column
#'
#' @return
#'
#' @usage
#' pam_from_table(data, ID_column, species_column)
#'
#' @export

pam_from_table <- function(data, ID_column, species_column) {
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





#' Creates a data frame of specie's references from RasterStack
#'
#' @description Creates a data frame of specie's references that contains longitude,
#' latitude, and species name, from a RasterStack or a RasterBrick.
#'
#' @param species_layers RasterStack or RasterBrick objects of species' presence (1)
#' absence (0)
#'
#' @return
#'
#' @usage
#' stack_2data(species_layers)
#'
#' @export
#' @importFrom raster rasterToPoints

stack_2data <- function(species_layers) {
  # Initial tests
  if (missing(species_layers)) {
    stop("Argument 'species_layers' must be defined")
  }

  # Stack to matrix
  sppm <- raster::rasterToPoints(species_layers)
  spnames <- colnames(sppm)[-c(1, 2)]

  # Preparing data
  sps <- lapply(1:length(spnames), function(x) {
    data.frame(sppm[sppm[, 2 + x] == 1, 1:2], spnames[x])
  })

  sps <- do.call(rbind, sps)
  colnames(sps) <- c("Longitude", "Latitude", "Species")

  return(sps)
}





#' Creates a data frame of specie's references from SpatialPolygonsDataFrame
#'
#' @description Creates a data frame of specie's references that contains identifiers
#' and species name, from a SpatialPolygonsDataFrame.
#'
#' @param spdf_object SpatialPolygonsDataFrame objects of species' presence (1)
#' absence (0)
#' @param spdf_grid Grids in the geography of the region of interest (output
#' of the function grid_from_region)
#'
#' @return
#'
#' @usage
#' spdf_2data(spdf_object, spdf_grid)
#'
#' @export
#' @importFrom sp over

spdf_2data <- function(spdf_object, spdf_grid) {
  # Initial tests
  if (missing(spdf_object)) {
    stop("Argument 'spdf_object' must be defined")
  }
  if (missing(spdf_grid)) {
    stop("Argument 'spdf_grid' must be defined")
  }

  # Names to be matched
  ID <- spdf_grid@data$ID
  spnames <- as.character(spdf_object@data$Species)

  # Preparing data
  sps <- lapply(1:length(spnames), function(x) {
    na.omit(data.frame(ID, sp::over(spdf_grid,
                                    spdf_object[spnames == spnames[x], ])))
  })

  sps <- do.call(rbind, sps)
  colnames(sps) <- c("ID", "Species")

  return(sps)
}





#' Creates a data frame of specie's references from a list
#'
#' @description Creates a data frame of specie's references that contains longitude,
#' latitude, and species name, from a list.
#'
#' @param raster_list list of RasterLayer objects of species' presence (1)
#' absence (0)
#'
#' @return
#'
#' @usage
#' rlist_2data(raster_list)
#'
#' @export
#' @importFrom raster raster

rlist_2data <- function(raster_list) {
  # Initial tests
  if (missing(raster_list)) {
    stop("Argument 'raster_list' must be defined")
  }

  # Running in loop for all elements of list
  sps <- lapply(1:length(spnames), function(x) {
    # raster to matrix
    sppm <- raster::rasterToPoints(raster_list[[x]])
    spname <- names(raster_list[[x]])

    # Preparing data
    data.frame(sppm[sppm[, 3] == 1, 1:2], spname)
  })

  sps <- do.call(rbind, sps)
  colnames(sps) <- c("Longitude", "Latitude", "Species")

  return(sps)
}




#'
#'
#' @description
#'
#' @param path
#' @param format
#' @param spdf_grid
#'
#' @return
#'
#' @usage
#' files_2data(path, format, spdf_grid = NULL)
#'
#' @export
#' @importFrom rgdal readOGR
#' @importFrom raster raster

files_2data <- function(path, format, spdf_grid = NULL) {
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

  # Finding files according to format
  if (format %in% c("shp", "gpkg")) {
    if (is.null(spdf_grid)) {
      stop("Argument 'spdf_grid' must be defined if 'format' is shp or gpkg")
    }
    # Names to be matched
    ID <- spdf_grid@data$ID

    if (format == "shp") {
      patt <- ".shp$"
      subs <- ".shp"
      mlist <- gsub(subs, "", list.files(path = path, pattern = patt))
      spnames <- mlist
    } else {
      patt <- ".gpkg$"
      subs <- ".gpkg"
      mlist <- list.files(path = path, pattern = patt)
      spnames <- gsub(subs, "", mlist)
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

  # Running in loop for all elements of list
  sps <- lapply(1:length(spnames), function(x) {
    if (format %in% c("shp", "gpkg")) {
      if (format == "shp") {
        rs <- rgdal::readOGR(dsn = path, layer = mlist[x], verbose = FALSE)
      } else {
        rs <- rgdal::readOGR(paste0(path, "/", mlist[x]), spnames[x],
                             verbose = FALSE)
      }

      sppm <- na.omit(data.frame(ID, Species = sp::over(spdf_grid, rs)[, 1]))
      sppm$Species <- spnames[x]
      return(sppm)

    } else {
      # Raster from file
      rs <- raster::raster(mlist[x])

      # Raster to matrix
      sppm <- raster::rasterToPoints(rs[[x]])

      # Preparing data
      return(data.frame(sppm[sppm[, 3] == 1, 1:2], spnames[x]))
    }
  })

  sps <- do.call(rbind, sps)

  if (format %in% c("shp", "gpkg")) {
    colnames(sps) <- c("ID", "Species")
  } else {
    colnames(sps) <- c("Longitude", "Latitude", "Species")
  }

  return(sps)
}
