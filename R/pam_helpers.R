#' Create the grids in the geography to be used to extract the values of the base_pam
#'
#' @description

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





stack_2data <- function(species_layers) {
  # Initial tests
  if (missing(species_layers)) {
    stop("Argument 'species_layers' must be defined")
  }

  # Stack to matrix
  sppm <- rasterToPoints(species_layers)
  spnames <- colnames(sppm)[-c(1,2)]

  # Preparing data
  sps <- lapply(1:length(spnames), function(x) {
    data.frame(sppm[sppm[, 2 + x] == 1, 1:2], spnames[x])
  })

  sps <- do.call(rbind, sps)
  colnames(sps) <- c("Longitude", "Latitude", "Species")

  return(sps)
}
