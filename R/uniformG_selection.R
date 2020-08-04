#' Selection of survey sites maximizing uniformity in geography
#'
#' @description Selection of sites to be sampled in a survey, with the goal of
#' maximizing uniformity of points in geographic space.
#'
#' @param master a master_matrix object derived from function
#' \code{\link{master_matrix}} or a master_selection object derived from functions
#' \code{\link{random_selection}}, \code{\link{uniformE_selection}},
#' or \code{\link{EG_selection}}.
#' @param expected_points (numeric) number of survey points (sites) to be selected.
#' @param max_n_samplings (numeric) maximum number of samples to be chosen after
#' performing all thinning \code{replicates}. Default = 1.
#' @param initial_distance (numeric) distance in km to be used for a first
#' process of thinning and detection of remaining points.
#' @param increase (numeric) value to be added to \code{initial_distance} until
#' reaching the number of \code{expected_points}.
#' @param replicates (numeric) number of thinning replicates. Default = 10.
#' @param use_preselected_sites (logical) whether to use sites that have been
#' defined as part of the selected sites previous any selection. Object in
#' \code{master} must contain the site(s) preselected in and element of name
#' "preselected_sites" for this argument to be effective. Default = TRUE.
#' See details for more information on the approach used.
#' @param n_optimization (numeric) number of times the algorithm of optimization
#' for site selection when using \code{use_preselected_sites} will run. Default
#' = 1000.
#' @param median_distance_filter (character) optional argument to define a median
#' distance-based filter based on which sets of sampling sites will be selected.
#' The default, NULL, does not apply such a filter. Options are: "max" and "min".
#' @param set_seed (numeric) integer value to specify a initial seed. Default = 1.
#'
#' @return
#' A master_selection object (S3) with an additional element called
#' selected_sites_G containing one or more sets of selected sites.
#'
#' @details
#' Survey sites are selected searching for maximum geographic distances among
#' all sites. This approach helps in selecting points that can cover most of the
#' geographic extent of the region of interest. This type of selection could be
#' appropriate when the region of interest has a complex geographic pattern (e.g.,
#' an archipelago). This type of selection does not consider environmental
#' conditions in the region of interest, which is why important environmental
#' combinations may not be represented in the final selection of sites.
#'
#' Exploring the geographic and environmental spaces of the region of interest
#' would be a crucial first step before selecting survey sites. Such explorations
#' can be done using the function \code{\link{explore_data_EG}}.
#'
#' If \code{use_preselected_sites} is TRUE and such sites are included as an
#' element in the object in \code{master}, the approach for selecting uniform
#' sites in geography is different than what was described above. User preselected
#' sites will always be part of the sites selected and other points are selected
#' based on an process of optimization that searches for the random set of points
#' that has the largest median distance among \code{n_optimization} distinct sets.
#'
#' As multiple sets could result from selection when the \code{use_preselected_sites}
#' is set as FALSE, the argument of the function \code{median_distance_filter}
#' could be used to select the set of sites with the maximum ("max") or minimum
#' ("min") median distance among selected sites. The option "max" will increase
#' the geographic distance among sampling sites, which could be desirable if the
#' goal is to cover the region of interest more broadly. The other option "min",
#' could be used in cases when the goal is to reduce resources and time needed
#' to sample such sites.
#'
#' @seealso
#' \code{\link{random_selection}}, \code{\link{uniformE_selection}},
#' \code{\link{EG_selection}}, \code{\link{plot_sites_EG}}
#'
#' @usage
#' uniformG_selection(master, expected_points, max_n_samplings = 1,
#'                    initial_distance, increase, replicates = 10,
#'                    use_preselected_sites = TRUE, n_optimization = 1000,
#'                    median_distance_filter = NULL, set_seed = 1)
#'
#' @export
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#'
#' # Selecting sites uniformly in G space
#' selectionG <- uniformG_selection(m_matrix, expected_points = 40, max_n_samplings = 1,
#'                                  initial_distance = 145, increase = 1,
#'                                  replicates = 5, set_seed = 1)

uniformG_selection <- function(master, expected_points, max_n_samplings = 1,
                               initial_distance, increase, replicates = 10,
                               use_preselected_sites = TRUE, n_optimization = 1000,
                               median_distance_filter = NULL, set_seed = 1) {
  # initial tests
  if (missing(master)) {
    stop("Argument 'master' is not defined.")
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
  if (!is.null(median_distance_filter)) {
    if (!median_distance_filter %in% c("max", "min")) {
      stop("Argument 'median_distance_filter' is not valid, see function's help.")
    }
  }
  if (use_preselected_sites == TRUE & is.null(master$preselected_sites)) {
    message("Element 'preselected_sites' in 'master' is NULL, setting\n'use_preselected_sites' = FALSE.")
    use_preselected_sites <- FALSE
  }

  # preparing data
  data <- master$master_matrix
  x_column <- "Longitude"
  y_column <- "Latitude"
  data <- data[!is.na(data[, x_column]) & !is.na(data[, y_column]), ]

  # selection depending on option of user points
  if (use_preselected_sites == TRUE) {
    # using preselected sites, optimization based on random selection
    n <- nrow(data)
    pre <- master$preselected_sites

    selected_sites <- lapply(1:n_optimization, function(x) {
      set.seed(set_seed + x - 1)
      sam <- sample(n, expected_points)
      dat <- unique(rbind(pre[, -1], data[sam, ]))
      dat[1:expected_points, ]
    })

    # Post filtering to get higher uniformity in G
    selected_sites <- distance_filter(selected_sites, "max")[1]

    # Returning results
    names(selected_sites) <- paste0("selection_", 1:length(selected_sites))
    master$selected_sites_G <- selected_sites

    return(structure(master, class = "master_selection"))

  } else {
    # preparing selection variables
    np <- nrow(data)
    ininp <- np
    dist <- initial_distance
    inin <- 1
    count <- 1

    # condition
    if (np < expected_points) {
      stop("Number of points in 'master_matrix' is smaller than 'expected_points'.")
    }
    if (np == expected_points) {
      message("Number of points in 'master_matrix' equals 'expected_points'.")
      master$selected_sites_G <- list(selection_1 = data)
      return(structure(master, class = "master_selection"))
    }

    # selection process
    while (np > expected_points) {
      # thinning
      thin <- point_thinning(data, x_column, y_column, dist, space = "G",
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
          master$selected_sites_G <- thin
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
}
