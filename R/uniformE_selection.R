#' Selection of survey sites maximizing uniformity in environmental space
#'
#' @description Selection of sites to be sampled in a survey, with the goal of
#' maximizing uniformity of points in environmental space.
#'
#' @param master master_matrix object derived from function
#' \code{\link{prepare_master_matrix}} or master_selection object derived
#' from functions \code{\link{random_selection}},
#' \code{\link{uniformG_selection}}, or \code{\link{EG_selection}}.
#' @param variable_1 (character or numeric) name or position of the first
#' variable (x-axis). If the function \code{\link{make_blocks}} was used in a
#' previous step, the default, NULL, will use the same two variables, otherwise
#' this argument must be defined.
#' @param variable_2 (character or numeric) name or position of the second
#' variable (y-axis). If the function \code{\link{make_blocks}} was used in a
#' previous step, the default, NULL, will use the same two variables, otherwise
#' this argument must be defined.
#' @param selection_from (character) set of points to perform the selection
#' from. Two options are available, "all_points" or "block_centroids". The
#' first option picks the points from all points in the environmental cloud,
#' and the second one selects points only from centroids of environmental
#' blocks. See \code{\link{make_blocks}}. Default = "all_points".
#' @param expected_points (numeric) total number of survey points (sites) to be
#' selected.
#' @param guess_distances (logical) whether or not to use internal algorithm
#' to automatically select \code{initial_distance} and \code{increase}. Default
#' = TRUE. If FALSE, \code{initial_distance} and \code{increase} must be
#' defined.
#' @param initial_distance (numeric) Euclidean distance to be used for a first
#' process of thinning and detection of remaining points. Default = NULL.
#' @param increase (numeric) initial value to be added to or subtracted from
#' \code{initial_distance} until reaching the number of \code{expected_points}.
#' Default = NULL.
#' @param max_n_samplings (numeric) maximum number of samples to be chosen
#' after performing all thinning \code{replicates}. Default = 1.
#' @param replicates (numeric) number of thinning replicates. Default = 10.
#' @param use_preselected_sites (logical) whether to use sites that have been
#' defined as part of the selected sites previous any selection. Object in
#' \code{master} must contain the site(s) preselected in and element of name
#' "preselected_sites" for this argument to be effective. Default = TRUE.
#' See details for more information on the approach used.
#' @param median_distance_filter (character) optional argument to define a
#' median distance-based filter based on which sets of sampling sites will be
#' selected. The default, NULL, does not apply such a filter. Options are:
#' "max" and "min".
#' @param set_seed (numeric) integer value to specify a initial seed.
#' Default = 1.
#' @param verbose (logical) whether or not to print messages about the process.
#' Default = TRUE.
#' @param force (logical) whether to replace existing set of sites selected
#' with this method in \code{master}.
#'
#' @return
#' A \code{\link{master_selection}} object (S3) with an element called
#' selected_sites_E containing one or more sets of selected sites.
#'
#' @details
#' Survey sites are selected in ways in which points will be uniformly
#' dispersed in environmental space, helping to select sites that present
#' different environmental conditions across the area of interest. This type of
#' selection is very useful to include, in the selected sites, distinct
#' environmental combinations existent in the area of interest. However, as the
#' distribution of climatic or other environmental combinations is not uniform
#' in geography, the sites selected with this function could appear clustered
#' when looked in a map.
#'
#' Exploring the geographic and environmental spaces of the region of interest
#' would be a crucial first step before selecting survey sites. Such
#' explorations can be done using the function \code{\link{explore_data_EG}}.
#'
#' If \code{use_preselected_sites} = TRUE and such sites are included as an
#' element in the object in \code{master}, the approach for selecting uniform
#' sites in environmental space is different than what was described above.
#' User-preselected sites will always be part of the sites selected. Other
#' points are selected based on an algorithm that searches for sites that are
#' uniformly distributed in environmental space but at a distance from
#' preselected sites that helps in maintaining uniformity. Note that
#' preselected sites will not be processed; therefore, uniformity of such points
#' cannot be warrantied. As multiple sets could result from selection, the
#' argument of the function \code{median_distance_filter} could be used to
#' select the set of sites with the maximum ("max") or minimum ("min") median
#' distance among selected sites. Option "max" will increase the geographic
#' distance among sampling sites, which could be desirable if the goal is to
#' cover the region of interest more broadly. The other option, "min", could be
#' used in cases when the goal is to reduce resources and time needed to sample
#' such sites.
#'
#' @seealso
#' \code{\link{uniformG_selection}}, \code{\link{random_selection}},
#' \code{\link{EG_selection}}, \code{\link{make_blocks}},
#' \code{\link{plot_sites_EG}}
#'
#' @usage
#' uniformE_selection(master, variable_1 = NULL, variable_2 = NULL,
#'                    selection_from = "all_points", expected_points,
#'                    guess_distances = TRUE, initial_distance = NULL,
#'                    increase = NULL, max_n_samplings = 1,
#'                    replicates = 10, use_preselected_sites = TRUE,
#'                    median_distance_filter = NULL, set_seed = 1,
#'                    verbose = TRUE, force = FALSE)
#'
#' @importFrom terra vect mask as.data.frame crs
#'
#' @export
#'
#' @examples
#' # Data
#' m_matrix <- read_master(system.file("extdata/m_matrix.rds",
#'                                     package = "biosurvey"))
#'
#' # Making blocks for analysis
#' m_blocks <- make_blocks(m_matrix, variable_1 = "PC1",
#'                         variable_2 = "PC2", n_cols = 10, n_rows = 10,
#'                         block_type = "equal_area")
#'
#' # Checking column names
#' colnames(m_blocks$data_matrix)
#'
#' # Selecting sites uniformly in E space
#' # because the make_blocks function was used, the same two variables will be
#' # used by default.
#' selectionE <- uniformE_selection(m_blocks, selection_from = "block_centroids",
#'                                  expected_points = 15, max_n_samplings = 1,
#'                                  replicates = 5, set_seed = 1)

uniformE_selection <- function(master, variable_1 = NULL, variable_2 = NULL,
                               selection_from = "all_points", expected_points,
                               guess_distances = TRUE, initial_distance = NULL,
                               increase = NULL, max_n_samplings = 1,
                               replicates = 10, use_preselected_sites = TRUE,
                               median_distance_filter = NULL, set_seed = 1,
                               verbose = TRUE, force = FALSE) {
  # Initial tests
  if (missing(master)) {
    stop("Argument 'master' is not defined.")
  }
  clsm <- class(master)[1]
  if (clsm %in% c("master_matrix", "master_selection")) {
    if (clsm == "master_selection") {
      if (!is.null(master$selected_sites_E) & force == FALSE) {
        stop("'master' already contains a selection of this type, use 'force' = TRUE to replace it")
      }
    }
  } else {
    stop("Argument 'master' must be of class 'master_matrix' or 'master_selection'")
  }
  if (is.null(master$data_matrix$Block)) {
    if (is.null(variable_1)) {
      stop("Argument 'variable_1' is not defined")
    }
    if (is.null(variable_2)) {
      stop("Argument 'variable_2' is not defined")
    }
  } else {
    sel_args <- attributes(master$data_matrix)

    variable_1 <- sel_args$arguments$variable_1
    variable_2 <- sel_args$arguments$variable_2
  }
  coln <- colnames(master$data_matrix)
  if (!variable_1 %in% coln) {
    stop(variable_1, " is not one o the columns in 'master$data_matrix'.")
  }
  if (!variable_2 %in% coln) {
    stop(variable_2, " is not one o the columns in 'master$data_matrix'.")
  }
  if (missing(expected_points)) {
    stop("Argument 'expected_points' is not defined.")
  }
  if (guess_distances == FALSE) {
    if (missing(initial_distance)) {
      stop("Argument 'initial_distance' is not defined.")
    }
    if (missing(increase)) {
      stop("Argument 'increase' is not defined.")
    }
  }
  if (max_n_samplings > replicates) {
    stop("Argument 'replicates' must be larger than 'max_n_samplings'.")
  }
  if (use_preselected_sites == TRUE & is.null(master$preselected_sites)) {
    if (verbose == TRUE) {
      message("Element 'preselected_sites' in 'master' is NULL, setting\n'use_preselected_sites' = FALSE")
    }
    use_preselected_sites <- FALSE
  }

  # Arguments for attributes
  other_args <- list(arguments = list(variable_1 = variable_1,
                                      variable_2 = variable_2,
                                      selection_from = selection_from,
                                      expected_points = expected_points,
                                      guess_distances = guess_distances,
                                      max_n_samplings = max_n_samplings,
                                      replicates = replicates,
                                      use_preselected_sites = use_preselected_sites,
                                      median_distance_filter = median_distance_filter,
                                      set_seed = set_seed))

  # Preparing data
  if (!selection_from[1] %in% c("all_points", "block_centroids")) {
    stop("Argument 'selection_from' is not valid, check function's help.")
  } else {
    # Preparing data
    data <- master$data_matrix

    if (selection_from[1] == "block_centroids") {
      if (is.null(master$data_matrix$Block)) {
        stop("Blocks are not defined in data_matrix, see function 'make_blocks'.")
      }

      # Preparing centroids
      data <- closest_to_centroid(data, variable_1, variable_2, space = "E",
                                  n = 1, id_column = "Block")
      useb <- TRUE
    } else {
      useb <- FALSE
    }
  }
  if (!is.null(median_distance_filter)) {
    if (!median_distance_filter %in% c("max", "min")) {
      stop("Argument 'median_distance_filter' is not valid, see function's help.")
    }
  }

  # Selection depending on option of user points
  if (use_preselected_sites == TRUE) {
    # Using preselected sites to create mask and define distance
    tst <- preselected_dist_mask(master, expected_points = expected_points,
                                 space = "E", variable_1 = variable_1,
                                 variable_2 = variable_2, use_blocks = useb,
                                 verbose = verbose)

    # Excluding close points from analysis
    data <- terra::vect(data, geom = c(x = variable_1, y = variable_2),
                        crs = terra::crs("+init=epsg:4326"))
    data <- terra::mask(data, tst$mask, inverse = TRUE)

    # Npre
    npre <- nrow(master$preselected_sites)
    expected_points <- expected_points - npre
  }

  # Preparing selection variables
  #Get data again in dataframe
  data <- terra::as.data.frame(data, geom = "XY")
  #Change xy to variables names
  colnames(data)[colnames(data) %in% c("x", "y")] <- c(variable_1, variable_2)

  ## Guessing distances
  if (guess_distances == TRUE) {
    if (use_preselected_sites == TRUE) {
      dist <- round(tst$distance, 2)
    } else {
      ext <- apply(data[, c(variable_1, variable_2)], 2, range)
      dist <- c(stats::dist(ext)) / expected_points
    }
    increase <- dist / 10
  } else {
    dist <- initial_distance
  }

  ## Other variables
  np <- nrow(data)
  ininp <- np
  pnp <- np
  inin <- 1
  count <- 1

  # Condition
  mess <- ifelse(selection_from[1] == "all_points",
                 "Number of available points using 'all_points'",
                 "Number of available points using 'block_centroids.
                 '")
  if (np < expected_points) {
    stop(mess, " is smaller than 'expected_points'.")
  }
  if (np == expected_points) {
    if (verbose == TRUE) {
      message(mess, " equals 'expected_points'")
    }

    if (use_preselected_sites == TRUE) {
      data <- rbind(master$preselected_sites[, -1],
                    data[, colnames(data) != "Selected_blocks"])
    }

    thin <- list(selection_1 = data)

    attributes(thin) <- c(attributes(thin), other_args)

    return(new_master_selection(master$data_matrix, master$preselected_sites,
                                master$region, master$mask, master$raster_base,
                                master$PCA_results, master$selected_sites_random,
                                master$selected_sites_G, thin,
                                master$selected_sites_EG))
  }

  # Selection process
  if (verbose == TRUE) {
    message("Running algorithm for selecting sites, please wait...")
  }

  while (np > expected_points) {
    # Thinning
    thin <- point_thinning(data, variable_1, variable_2, dist, space = "E",
                           max_n_samplings, replicates, set_seed)
    np <- nrow(thin[[1]])
    if (verbose == TRUE) {
      message("    Distance  ", round(dist, 3), "  resulted in  ",
              np, "  points")
    }

    if (np <= expected_points) {
      if (np == expected_points) {
        # Selection done
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

  # Returning results
  ## Adding preselected to selected
  if (use_preselected_sites == TRUE) {
    thin <- lapply(thin, function(x) {
      rbind(master$preselected_sites[, -1], x[, colnames(x) != "Selected_blocks"])
    })
  }

  ## Filtering by distances success
  if (length(thin) > 1 & !is.null(median_distance_filter)) {
    thin <- distance_filter(thin, median_distance_filter)
  }

  ## Naming and returning
  names(thin) <- paste0("selection_", 1:length(thin))

  if (verbose == TRUE) {
    message("Total number of sites selected: ", nrow(thin[[1]]))
  }

  attributes(thin) <- c(attributes(thin), other_args)

  if (class(master)[1] == "master_matrix") {
    return(new_master_selection(master$data_matrix, master$preselected_sites,
                                master$region, master$mask, master$raster_base,
                                master$PCA_results, selected_sites_random = NULL,
                                selected_sites_G = NULL, thin,
                                selected_sites_EG = NULL))

  } else {
    return(new_master_selection(master$data_matrix, master$preselected_sites,
                                master$region, master$mask, master$raster_base,
                                master$PCA_results, master$selected_sites_random,
                                master$selected_sites_G, thin,
                                master$selected_sites_EG))
  }
}
