#' Selection of survey sites maximizing uniformity in environmental space
#'
#' @description Selection of sites to be sampled in a survey, with the goal of
#' maximizing uniformity of points in environmental space.
#'
#' @param master a master_matrix object derived from function
#' \code{\link{master_matrix}} or a master_selection object derived from functions
#' \code{\link{random_selection}}, \code{\link{uniformG_selection}},
#' or \code{\link{EG_selection}}.
#' @param variable_1 (character or numeric) name or position of the first
#' variable (X-axis).
#' @param variable_2 (character or numeric) name or position of the second
#' variable (Y-axis).
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
#' A master_selection object (S3) with an aditional element called
#' selected_sites_E containing one or more sets of selected sites.
#'
#' @usage
#' uniformE_selection(master, variable_1, variable_2, selection_from = "all_points",
#'                    expected_points, max_n_samples = 1,
#'                    initial_distance, increase, replicates = 10, set_seed = 1)
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
#' selectionE <- uniformE_selection(m_matrix, variable_1 = "PC1", variable_2 = "PC2",
#'                                  selection_from = "block_centroids",
#'                                  expected_points = 15, max_n_samples = 1,
#'                                  initial_distance = 1, increase = 0.1,
#'                                  replicates = 5, set_seed = 1)

uniformE_selection <- function(master, variable_1, variable_2,
                               selection_from = "all_points", expected_points,
                               max_n_samples = 1, initial_distance, increase,
                               replicates = 10, set_seed = 1) {
  # Initial tests
  if (missing(master)) {
    stop("Argument 'master' is not defined.")
  }
  if (missing(variable_1)) {
    stop("Argument 'variable_1' is not defined.")
  }
  if (missing(variable_2)) {
    stop("Argument 'variable_2' is not defined.")
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
      data <- master$master_matrix
      data <- data[!is.na(data[, variable_1]) & !is.na(data[, variable_2]), ]
    } else {
      data <- master$master_matrix

      # preparing centroids
      data <- closest_to_centroid(data, variable_1, variable_2, space = "E",
                                  n = 1, id_column = "Block")
    }
  }

  # preaparing selection variables
  np <- nrow(data)
  ininp <- np
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
    master$selected_sites_E <- list(selection_1 = data)
    return(structure(master, class = "master_selection"))
  }

  # selection process
  while (np > expected_points) {
    # thinning
    thin <- point_thinning(data, variable_1, variable_2, dist, space = "E",
                           max_n_samples, replicates, set_seed)
    np <- nrow(thin[[1]])
    message("    Distance  ", dist, "  resulted in  ", np, "  points")

    if (np <= expected_points) {
      if (np == expected_points) {
        # success
        names(thin) <- paste0("selection_", 1:length(thin))
        master$selected_sites_E <- thin
        return(structure(master, class = "master_selection"))
      } else {
        if (count == 1) {
          stop("'initial_distance' resulted in  ", np, "  points. Try smaller values.")
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

          np <- ininp
        }
      }
    }

    # starting again
    count <- count + 1
    pdist <- dist
    dist <- dist + increase
  }
}
