#' Create the grids in the geography to be used to extract the values of the base_pam
#'
#' @description
#'
#' @param region SpatialPolygonsDataFrame of the region of interest.
#' @param cell_size
#'
#' @return
#'
#' @usage
#' grid_from_region(region, cell_size)
#'
#' @export
#' @importFrom raster extent raster res values mask rasterToPolygons rasterToPoints
#' @importFrom sp proj4string
#'
#' @examples

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





#'
#'
#' @description
#'
#' @param data
#' @param ID_column
#' @param species_column
#'
#' @return
#'
#' @usage
#' pam_from_table(data, ID_column, species_column)
#'
#' @export
#'
#' @examples

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
  counts <- data.frame(rownames(counts), counts) # as.numeric(rownames(counts)
  colnames(counts) <- c(ID_column, nams)

  return(counts)
}





#'
#'
#' @description
#'
#' @param species_layers
#'
#' @return
#'
#' @usage
#' stack_2data(species_layers)
#'
#' @export
#' @importFrom raster rasterToPoints
#'
#' @examples

stack_2data <- function(species_layers) {
  # Initial tests
  if (missing(species_layers)) {
    stop("Argument 'species_layers' must be defined")
  }

  # Stack to matrix
  sppm <- raster::rasterToPoints(species_layers)
  spnames <- colnames(sppm)[-c(1,2)]

  # Preparing data
  sps <- lapply(1:length(spnames), function(x) {
    data.frame(sppm[sppm[, 2 + x] == 1, 1:2], spnames[x])
  })

  sps <- do.call(rbind, sps)
  colnames(sps) <- c("Longitude", "Latitude", "Species")

  return(sps)
}





#'
#'
#' @description
#'
#' @param species_layers
#'
#' @return
#'
#' @usage
#' spdf_2data(spdf_object, spdf_grid)
#'
#' @export
#' @importFrom sp over
#'
#' @examples

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
