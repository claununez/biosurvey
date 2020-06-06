#' Selection of survey sites maximizing uniformity in environmental space
#'
#' @description Selection of sites to be sampled in a survey, with the goal of
#' maximizing uniformity of points in environmental space.
#'
#' @param master_matrix object derived from function \code{\link{master_matrix}}.
#' Optionally, if master_matrix is not necessary, a list containing an object of
#' class data.frame with at least two columns represening two variables. The name
#' of this element in the list must be "master_matrix". For instance:
#' \code{my_list <- list(master_matrix = YOUR_data.frame)}.
#' @param x_column (character) the name of the X-axis.
#' @param y_column (character) the name of the Y-axis.
#' @param selection_from (character) set of points to perfomr the selection from.
#' Two options are available, "all_points" or "block_centroids". The first option
#' picks the points from all points in the environmental cloud, and the second
#' one selects points only from centroids of environmental blocks. See
#' \code{\link{make_blocks}}. Default = "all_points".
#' @param expected_points (numeric) number of survey points (sites) to be selected.
#' @param max_n_samples (numeric) maximun number of samples to be chosen after
#' performing all thinning \code{replicates}. Default = 1.
#' @param initial_distance (numeric) euclidean distance to be used for a first
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
#' uniformE_selection(master_matrix, x_column, y_column,
#'                    selection_from = "all_points", expected_points,
#'                    max_n_samples = 1, initial_distance, increase,
#'                    replicates = 10, set_seed = 1)
#'
#' @export
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#'
#' # Making blocks for analysis
#' m_matrix <- make_blocks(m_matrix, variable_1 = "PC1",
#'                         variable_2 = "PC2", n_cols = 10, n_rows = 10,
#'                         block_type = "equal_area")
#'
#' # Checking column names
#' colnames(m_matrix$master_matrix)
#' summary(m_matrix$master_matrix[, 9:10])
#'
#' # Selecting sites uniformly in E space
#' selectionE <- uniformE_selection(m_matrix, x_column = "PC1", y_column = "PC2",
#'                                  selection_from = "block_centroids",
#'                                  expected_points = 15, max_n_samples = 1,
#'                                  initial_distance = 1, increase = 0.1,
#'                                  replicates = 5, set_seed = 1)

uniformE_selection <- function(master_matrix, x_column, y_column,
                               selection_from = "all_points", expected_points,
                               max_n_samples = 1, initial_distance, increase,
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
  if (!selection_from[1] %in% c("all_points", "block_centroids")) {
    stop("Argument 'selection_from' is not valid, check function's help.")
  } else {
    # preparing data
    if (selection_from[1] == "all_points") {
      data <- master_matrix$master_matrix
      data <- data[!is.na(data[, x_column]) & !is.na(data[, y_column]), ]
    } else {
      data <- master_matrix$master_matrix

      # preparing centroids
      data <- closest_to_centroid(data, x_column, y_column, space = "E",
                                  n = 1, id_column = "Block")
    }
  }

  # preaparing selection variables
  np <- nrow(data)
  dist <- initial_distance
  inin <- 1
  count <- 1

  # condition
   mess <- ifelse(selection_from[1] == "all_points",
                 "Number of points in 'master_matrix'",
                 "Number of block centroid points")
  if (np < expected_points) {
    stop(mess, " is smaller than 'expected_points'.")
  }
  if (np == expected_points) {
    message(mess, " equals 'expected_points'.")
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
        names(thin) <- paste0("selection_", 1:length(thin))
        master_matrix$selected_sites <- thin
        return(master_matrix)
      } else {
        if (count == 1) {
          stop("'initial_distance' resulted in  ", np, "  points. Try smaler values.")
        } else {
          # reducing increase distance
          message("\nNo distance resulted in  ", expected_points, "  points.",
                  "\nTrying distances between  ", pdist, "  and  ", dist, ".\n")
          dist <- pdist

          increase <- increase / 10
          inin <- inin + 1

          if (inin == 3) {
            stop("\nNo distance resulted in  ", expected_points,
                 "  points after trying smaller intervals.",
                 "\nTry a distance of  ", pdist,
                 "  with values for 'increase' below  ", increase, ".")
          }
        }
      }
    }

    # starting again
    count <- count + 1
    pdist <- dist
    dist <- dist + increase
  }
}
