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
#' @param selection_from (character) set of points to perform the selection from.
#' Two options are available, "all_points" or "block_centroids". The first option
#' picks the points from all points in the environmental cloud, and the second
#' one selects points only from centroids of environmental blocks. See
#' \code{\link{make_blocks}}. Default = "all_points".
#' @param expected_points (numeric) number of survey points (sites) to be selected.
#' @param max_n_samplings (numeric) maximum number of samples to be chosen after
#' performing all thinning \code{replicates}. Default = 1.
#' @param initial_distance (numeric) euclidean distance to be used for a first
#' process of thinning and detection of remaining points.
#' @param increase (numeric) value to be added to \code{initial_distance} until
#' reaching the number of \code{expected_points}.
#' @param replicates (numeric) number of thinning replicates. Default = 10.
#' @param median_distance_filter (character) optional argument to define a median
#' distance-based filter based on which sets of sampling sites will be selected.
#' The default, NULL, does not apply such a filter. Options are: "max" and "min".
#' @param set_seed (numeric) integer value to specify a initial seed. Default = 1.
#'
#' @return
#' A master_selection object (S3) with an additional element called
#' selected_sites_E containing one or more sets of selected sites.
#'
#' @details
#' Survey sites are selected in ways in which points will be uniformly dispersed
#' in environmental space, helping to select sites that present different
#' environmental conditions across the area of interest. This type of selection
#' is very useful to include, in the selected sites, distinct environmental
#' combinations existent in the area of interest. However, as the distribution of
#' climatic or other environmental combinations is not uniform in geography, the
#' sites selected with this function could appear clustered when looked in a map.
#'
#' Exploring the geographic and environmental spaces of the region of interest
#' would be a crucial first step before selecting survey sites. Such explorations
#' can be done using the function \code{\link{explore_data_EG}}.
#'
#' @seealso
#' \code{\link{uniformG_selection}}, \code{\link{random_selection}},
#' \code{\link{EG_selection}}, \code{\link{make_blocks}},
#' \code{\link{plot_sites_EG}}
#'
#' @usage
#' uniformE_selection(master, variable_1, variable_2, selection_from = "all_points",
#'                    expected_points, max_n_samplings = 1,
#'                    initial_distance, increase, replicates = 10,
#'                    median_distance_filter = NULL, set_seed = 1)
#'
#' @export
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#'
#' # Making blocks for analysis
#' m_blocks <- make_blocks(m_matrix, variable_1 = "PC1",
#'                         variable_2 = "PC2", n_cols = 10, n_rows = 10,
#'                         block_type = "equal_area")
#'
#' # Checking column names
#' colnames(m_blocks$master_matrix)
#'
#' # Selecting sites uniformly in E space
#' selectionE <- uniformE_selection(m_blocks, variable_1 = "PC1", variable_2 = "PC2",
#'                                  selection_from = "block_centroids",
#'                                  expected_points = 15, max_n_samplings = 1,
#'                                  initial_distance = 1, increase = 0.1,
#'                                  replicates = 5, set_seed = 1)

uniformE_selection <- function(master, variable_1, variable_2,
                               selection_from = "all_points", expected_points,
                               max_n_samplings = 1, initial_distance, increase,
                               replicates = 10, median_distance_filter = NULL,
                               set_seed = 1) {
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
  if (!is.null(median_distance_filter)) {
    if (!median_distance_filter %in% c("max", "min")) {
      stop("Argument 'median_distance_filter' is not valid, see function's help.")
    }
  }

  # preparing selection variables
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
                           max_n_samplings, replicates, set_seed)
    np <- nrow(thin[[1]])
    message("    Distance  ", dist, "  resulted in  ", np, "  points")

    if (np <= expected_points) {
      if (np == expected_points) {
        # success
        if (length(thin) > 1 & !is.null(median_distance_filter)) {
          thin <- distance_filter(thin, median_distance_filter)

        }
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
