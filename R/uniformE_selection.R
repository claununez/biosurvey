
uniformE_selection <- function(master_matrix, x_column, y_column,
                               selection_from = "all_points",
                               initial_distance, increase,
                               expected_points, max_n_samples = 1,
                               replicates = 10, set_seed = 1) {
  # Initial tests
  if (missing(master_matrix)) {
    stop("Argument 'master_matrix' is not defined.")
  }
  if (missing(x_column)) {
    stop("Argument 'x_column' is not defined.")
  }
  if (missing(y_column)) {
    stop("Argument 'y_column' is not defined.")
  }
  if (missing(initial_distance)) {
    stop("Argument 'initial_distance' is not defined.")
  }
  if (missing(increase)) {
    stop("Argument 'increase' is not defined.")
  }
  if (missing(expected_points)) {
    stop("Argument 'expected_points' is not defined.")
  }
  if (missing(space)) {
    stop("Argument 'space' is not defined.")
  }
  if (!selection_from[1] %in% c("all_points", "centroids")) {
    stop("Argument 'selection_from' is not valid, check function's help.")
  } else {
    # preparing data
    if (selection_from[1] == "all_points") {
      data <- master_matrix$master_matrix
      data <- data[!is.na(data[, x_column]) & !is.na(data[, y_column]), ]
    } else {
      data <- master_matrix$master_matrix

      # preparing centroids

    }
  }



  # preaparing selection variables
  np <- nrow(data)
  dist <- initial_distance
  inin <- 1

  # condition
  if (np < expected_points) {
    stop("Number of points in 'master_matrix' is smaller than 'expected_points'.")
  }
  if (np == expected_points) {
    message("Number of points in 'master_matrix' equals 'expected_points'.")
    master_matrix$selected_sites <- list(selection_1 = data)
    return(master_matrix)
  }

  # selection process
  while (np > expected_points) {
    # thinning
    thin <- point_thinning(data, x_column, y_column, dist, space = "E",
                           max_n_samples, replicates, set_seed)
    np <- nrow(thin[[1]])
    message("    Distance  ", dist, "  resulted in  ", np, "  points")

    if (np <= expected_points) {
      if (np == expected_points) {
        # success
        master_matrix$selected_sites <- thin
        return(master_matrix)
      } else {
        # reducing increase distance
        message(paste("No distance resulted in ", expected_points, " points.",
                      "\nTrying distances between:\t", pdist, "and", dist))
        dist <- pdist
        if (inin == 1) {
          increase <- increase / 10
          inin <- 2
        }
        if (inin == 2) {
          increase <- increase / 10
          inin <- 3
        }
        if (inin == 3) {
          stop(paste("No distance resulted in", expected_points,
                     "points after trying smaller intervals.",
                     "n\Try a distance of", pdist,
                     "with a values for 'increase' below", increase))
        }
      }
    }

    # starting again
    pdist <- dist
    dist <- dist + increase
  }
}
