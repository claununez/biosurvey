#' Selection of survey sites maximizing uniformity in geography
#'
#' @description Selection of sites to be sampled in a survey, with the goal of
#' maximizing uniformity of points in geographic space.
#'
#' @param master_matrix object derived from function \code{\link{master_matrix}}.
#' Optionally, if master_matrix is not necessary, a list containing an object of
#' class data.frame with at least two columns represening two variables. The name
#' of this element in the list must be "master_matrix". For instance:
#' \code{my_list <- list(master_matrix = YOUR_data.frame)}.
#' @param expected_points (numeric) number of survey points (sites) to be selected.
#' @param max_n_samples (numeric) maximun number of samples to be chosen after
#' performing all thinning \code{replicates}. Default = 1.
#' @param initial_distance (numeric) distance in km to be used as for a first
#' process of thinning and detection of remaining points.
#' @param increase (numeric) value to be added to \code{initial_distance} untill
#' reaching the number of \code{expected_points}.
#' @param replicates (numeric) number of thinning replicates. Default = 10.
#' @param set_seed (numeric) integer value to specify a initial seed. Default = 1.
#'
#' @return
#' The master_matrix list with an aditional element containing one or more sets
#' of selected sites.
#'
#' @usage
#' uniformG_selection(master_matrix, expected_points, max_n_samples = 1,
#'                    initial_distance, increase, replicates = 10, set_seed = 1)
#'
#' @export



uniformG_selection <- function(master_matrix, expected_points, max_n_samples = 1,
                               initial_distance, increase, replicates = 10,
                               set_seed = 1) {
  # initial tests
  if (missing(master_matrix)) {
    stop("Argument 'master_matrix' is not defined.")
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

  # preparing data
  data <- master_matrix$master_matrix
  x_column <- "Longitude"
  y_column <- "Latitude"
  data <- data[!is.na(data[, x_column]) & !is.na(data[, y_column]), ]

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
    thin <- point_thinning(data, x_column, y_column, dist, space = "G",
                           max_n_samples, replicates, set_seed)
    np <- nrow(thin[[1]])
    message("    Distance  ", dist, "  resulted in  ", np, "  points")

    if (np <= expected_points) {
      if (np == expected_points) {
        # success
        names(thin) <- paste0("selection_", 1:length(thin))
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
                     "\nTry a distance of", pdist,
                     "with a values for 'increase' below", increase))
        }
      }
    }
    # starting again
    pdist <- dist
    dist <- dist + increase
  }
}