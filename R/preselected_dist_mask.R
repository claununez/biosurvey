#' Helper to create objects to detect points to close to preselected sites
#'
#' @param master master_matrix object derived from function
#' \code{\link{prepare_master_matrix}} or master_selection object derived
#' from functions \code{\link{random_selection}},
#' \code{\link{uniformE_selection}}, or \code{\link{EG_selection}}.
#' @param expected_points (numeric) number of survey points (sites) to be
#' selected.
#' @param space (character) space in which the thinning will be performed. There
#' are two options available: "G", if it will be in geographic space, and
#' "E", if it will be in environmental space.
#' @param variable_1 (character or numeric) name or position of the first
#' variable (x-axis). Default = NULL.
#' @param variable_2 (character or numeric) name or position of the second
#' variable (y-axis). Default = NULL.
#' @param use_blocks (logical) whether or not to use block centroids instead of
#' all points when \code{space} = "E". Default = FALSE.
#' @param verbose (logical) whether or not to print messages about the process.
#' Default = TRUE.
#'
#' @return
#' A list of two elements: the distance used to obtain \code{expected_points}
#' and a SpatVector object created from preselected_sites in
#' master.
#'
#' @usage
#' preselected_dist_mask(master, expected_points, space, variable_1 = NULL,
#'                       variable_2 = NULL, use_blocks = FALSE,
#'                       verbose = TRUE)
#'
#' @export
#' @importFrom terra distance buffer disagg project crs vect
#' @importFrom stats dist
#'
#' @examples
#' # Data
#' data("m_matrix_pre", package = "biosurvey")
#'
#' # Running
#' pdm <- preselected_dist_mask(master = m_matrix_pre, expected_points = 20,
#'                              space = "G")

preselected_dist_mask <- function(master, expected_points, space,
                                  variable_1 = NULL, variable_2 = NULL,
                                  use_blocks = FALSE, verbose = TRUE) {

  # Initial tests
  if (missing(master)) {
    stop("Argument 'master' is not defined.")
  }
  if (missing(expected_points)) {
    stop("Argument 'expected_points' is not defined.")
  }
  if (missing(space)) {
    stop("Argument 'space' is not defined.")
  }
  if (space == "G" & use_blocks == TRUE) {
    use_blocks <- FALSE
  }
  if (space == "E" & any(is.null(c(variable_1, variable_2)))) {
    stop("Arguments 'variable_1' and 'variable_2' must be defined.")
  }

  # Preparing data
  data <- master$data_matrix
  pre <- master$preselected_sites

  if (use_blocks == TRUE) {
    if (is.null(master$data_matrix$Block)) {
      stop("Blocks are not defined in data_matrix, see function 'make_blocks'.")
    }

    # Preparing centroids
    data <- closest_to_centroid(data, variable_1, variable_2, space = "E",
                                n = 1, id_column = "Block")
  }

  # Test if enough points
  np <- nrow(data)

  ## Condition
  mess <- ifelse(use_blocks == FALSE, "Number of points in 'data_matrix'",
                 "Number of block centroid points")
  if (np < expected_points) {
    stop(mess, " is smaller than 'expected_points'.")
  }
  if (np == expected_points) {
    if (verbose == TRUE) {
      message("Number of points in 'data_matrix' equals 'expected_points'.")
    }
  }

  # Processing
  if (space == "G") {
    x_column <- "Longitude"
    y_column <- "Latitude"
    gv <- c(x_column, y_column)
    ext <- apply(data[, gv], 2, range)
    mxdis <- terra::distance(x = ext, lonlat = TRUE)
    dist <- (mxdis / 1000) / expected_points
    increase <- dist / 10
  } else {
    x_column <- variable_1
    y_column <- variable_2
    gv <- c(x_column, y_column)

    ext <- apply(data[, gv], 2, range)
    dist <- c(stats::dist(ext)) / expected_points
    increase <- dist / 10
  }

  # Preparing selection variables
  ininp <- np
  pnp <- np
  inin <- 1
  count <- 1

  # Testing distance
  if (verbose == TRUE) {
    message("Running distance optimization, please wait...")
  }

  while (np > expected_points) {
    # Thinning
    thin <- point_thinning(data, x_column, y_column, dist, space, 1, 1)
    np <- nrow(thin[[1]])

    if (np <= expected_points) {
      if (np == expected_points) {
        # Distance found
        break()
      } else {
        # Reducing initial distance
        dist <- dist - increase

        # Reducing increase distance
        if (count > 1 & pnp > expected_points) {
          increase <- increase / 10
        }
      }
    } else {
      dist <- dist + increase
    }

    # Starting again
    pnp <- np
    np <- ininp
    count <- count + 1
  }

  # Distance determined, creating mask
  if (space == "G") {
    sppre <- wgs84_2aed_laea(pre, x_column, y_column)
    maskp <- terra::buffer(sppre, width = dist * 1000, quadsegs = 100)
    maskp <- terra::disagg(maskp)
    maskp$ID <- 1:length(maskp)
    maskp <- terra::project(maskp, terra::crs("+init=epsg:4326"))
  } else {
    if (use_blocks == TRUE) {
      pre <- data.frame(data[data$Block %in% pre$Block, ])
    }
    sppre <- terra::vect(pre, geom = c(x = variable_1, y = variable_2),
                         crs = terra::crs("+init=epsg:4326"))
    maskp <- suppressWarnings(terra::buffer(sppre, width = dist * 1000,
                                            quadsegs = 100))
    maskp <- terra::disagg(maskp)
    maskp$ID <- 1:length(maskp)
  }

  # Returning results
  return(list(distance = dist, mask = maskp))
}

