preselected_dist_mask <- function(master, expected_points, space, variable_1 = NULL,
                                  variable_2 = NULL, use_blocks = FALSE) {

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
    # preparing centroids
    data <- closest_to_centroid(data, variable_1, variable_2, space = "E",
                                n = 1, id_column = "Block")
  }

  # Test if enough points
  np <- nrow(data)

  ## condition
  mess <- ifelse(use_blocks == FALSE, "Number of points in 'data_matrix'",
                 "Number of block centroid points")
  if (np < expected_points) {
    stop(mess, " is smaller than 'expected_points'.")
  }
  if (np == expected_points) {
    message("Number of points in 'data_matrix' equals 'expected_points'.")
  }

  # Processing
  if (space == "G") {
    x_column <- "Longitude"
    y_column <- "Latitude"
    gv <- c(x_column, y_column)

    alld <- apply(raster::pointDistance(data[, gv], pre[, gv], lonlat = TRUE),
                  1, min)

    ext <- apply(data[, gv], 2, range)
    mxdis <- raster::pointDistance(ext[1, ], ext[2, ], lonlat = TRUE)
    dist <- (mxdis / 1000) / expected_points
    increase <- dist / 10
  } else {
    x_column <- variable_1
    y_column <- variable_2
    gv <- c(x_column, y_column)

    ep <- t(data[, gv])
    pres <- t(pre[, gv])

    alld <- apply(apply(pres, 2, function(pres) {sqrt(colSums((ep - pres)^2))}),
                  1, min)

    ext <- apply(data[, gv], 2, range)
    dist <- c(stats::dist(ext)) / expected_points
    increase <- dist / 10
  }

  # preparing selection variables
  ininp <- np
  inin <- 1
  count <- 1

  # testing distance
  message("Running distance optimization, please wait...")

  while (np > expected_points) {
    # thinning
    thin <- point_thinning(data, x_column, y_column, dist, space, 1, 1)
    np <- nrow(thin[[1]])

    if (np <= expected_points) {
      if (np == expected_points) {
        # distance found
        break()
      } else {
        # reducing initial distance
        dist <- dist - increase

        # reducing increase distance
        if (count > 1 & pnp > expected_points) {
          increase <- increase / 10
        }
      }
    } else {
      dist <- dist + increase
    }

    # starting again
    pnp <- np
    np <- ininp
    count <- count + 1
  }

  # distance determined, creating mask
  sppre <- wgs84_2aed_laea(pre, x_column, y_column)

  maskp <- rgeos::gBuffer(sppre, width = dist * 1000, quadsegs = 100)
  maskp <- sp::spTransform(maskp, sp::CRS("+init=epsg:4326"))

  # returning results
  return(list(distance = dist, mask = maskp))
}
