#' Selection of survey sites maximizing uniformity in environmental space
#'
#' @description Selection of sites to be sampled in a survey, with the goal of
#' maximizing uniformity of points in environmental space.
#'
#' @param master a master_matrix object derived from function
#' \code{\link{prepare_master_matrix}} or a master_selection object derived from functions
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
#' @param guess_distances (logical) whether or not to use internal algorithm
#' to automatically select \code{initial_distance} and \code{increase}. Default
#' = TRUE. If FALSE, \code{initial_distance} and \code{increase} must be defined.
#' @param initial_distance (numeric) euclidean distance to be used for a first
#' process of thinning and detection of remaining points. Default = NULL.
#' @param increase (numeric) initial value to be added to or subtracted from
#' \code{initial_distance} until reaching the number of \code{expected_points}.
#' Default = NULL.
#' @param max_n_samplings (numeric) maximum number of samples to be chosen after
#' performing all thinning \code{replicates}. Default = 1.
#' @param replicates (numeric) number of thinning replicates. Default = 10.
#' @param use_preselected_sites (logical) whether to use sites that have been
#' defined as part of the selected sites previous any selection. Object in
#' \code{master} must contain the site(s) preselected in and element of name
#' "preselected_sites" for this argument to be effective. Default = TRUE.
#' See details for more information on the approach used.
#' @param median_distance_filter (character) optional argument to define a median
#' distance-based filter based on which sets of sampling sites will be selected.
#' The default, NULL, does not apply such a filter. Options are: "max" and "min".
#' @param set_seed (numeric) integer value to specify a initial seed. Default = 1.
#' @param verbose (logical) whether or not to print messages about the process.
#' Default = TRUE.
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
#' If \code{use_preselected_sites} is TRUE and such sites are included as an
#' element in the object in \code{master}, the approach for selecting uniform
#' sites in environmental space is different than what was described above.
#' User preselected sites will always be part of the sites selected. Other points
#' are selected based on an algorithm that searches for sites that are uniformly
#' distributed in environmental space but at a distance from preselected sites
#' that helps in maintaining uniformity. Note that preselected sites will not be
#' processed, therefore, uniformity of such points cannot be warrantied.
#'
#' As multiple sets could result from selection, the argument of the function
#' \code{median_distance_filter} could be used to select the set of sites with
#' the maximum ("max") or minimum ("min") median distance among selected sites.
#' Option "max" will increase the geographic distance among sampling sites, which
#' could be desirable if the goal is to cover the region of interest more broadly.
#' The other option "min", could be used in cases when the goal is to reduce
#' resources and time needed to sample such sites.
#'
#' @seealso
#' \code{\link{uniformG_selection}}, \code{\link{random_selection}},
#' \code{\link{EG_selection}}, \code{\link{make_blocks}},
#' \code{\link{plot_sites_EG}}
#'
#' @usage
#' uniformE_selection(master, variable_1, variable_2,
#'                    selection_from = "all_points", expected_points,
#'                    guess_distances = TRUE, initial_distance = NULL,
#'                    increase = NULL, max_n_samplings = 1,
#'                    replicates = 10, use_preselected_sites = TRUE,
#'                    median_distance_filter = NULL, set_seed = 1,
#'                    verbose = TRUE)
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
#' colnames(m_blocks$data_matrix)
#'
#' # Selecting sites uniformly in E space
#' selectionE <- uniformE_selection(m_blocks, variable_1 = "PC1", variable_2 = "PC2",
#'                                  selection_from = "block_centroids",
#'                                  expected_points = 15, max_n_samplings = 1,
#'                                  initial_distance = 1, increase = 0.1,
#'                                  replicates = 5, set_seed = 1)

uniformE_selection <- function(master, variable_1, variable_2,
                               selection_from = "all_points", expected_points,
                               guess_distances = TRUE, initial_distance = NULL,
                               increase = NULL, max_n_samplings = 1,
                               replicates = 10, use_preselected_sites = TRUE,
                               median_distance_filter = NULL, set_seed = 1,
                               verbose = TRUE) {
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
      message("Element 'preselected_sites' in 'master' is NULL, setting\n'use_preselected_sites' = FALSE.")
    }
    use_preselected_sites <- FALSE
  }

  # preparing data
  if (!selection_from[1] %in% c("all_points", "block_centroids")) {
    stop("Argument 'selection_from' is not valid, check function's help.")
  } else {
    # preparing data
    data <- master$data_matrix

    if (selection_from[1] == "block_centroids") {
      if (is.null(master$data_matrix$Block)) {
        stop("Blocks are not defined in data_matrix, see function 'make_blocks'.")
      }

      # preparing centroids
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

  # selection depending on option of user points
  if (use_preselected_sites == TRUE) {
    # using preselected sites to create mask and define distance
    tst <- preselected_dist_mask(master, expected_points = expected_points,
                                 space = "E", variable_1 = variable_1,
                                 variable_2 = variable_2, use_blocks = useb,
                                 verbose = verbose)

    # excluding close points from analysis
    data <- sp::SpatialPointsDataFrame(data[, c(variable_1, variable_2)], data,
                                       proj4string = sp::CRS("+init=epsg:4326"))
    data <- data[is.na(sp::over(data, tst$mask)$ID), ]@data

    # npre
    npre <- nrow(master$preselected_sites)
    expected_points <- expected_points - npre
  }

  # preparing selection variables
  ## guessing distances
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

  ## other variables
  np <- nrow(data)
  ininp <- np
  pnp <- np
  inin <- 1
  count <- 1

  # condition
  mess <- ifelse(selection_from[1] == "all_points",
                 "Number of available points using 'all_points'",
                 "Number of available points using 'block_centroids.
                 '")
  if (np < expected_points) {
    stop(mess, " is smaller than 'expected_points'.")
  }
  if (np == expected_points) {
    if (verbose == TRUE) {
      message(mess, " equals 'expected_points'.")
    }

    if (use_preselected_sites == TRUE) {
      data <- rbind(master$preselected_sites[, -1],
                    data[, colnames(data) != "Selected_blocks"])
    }

    return(new_master_selection(master$data_matrix, master$preselected_sites,
                                master$region, master$mask, master$raster_base,
                                master$PCA_results, master$selected_sites_random,
                                master$selected_sites_G, list(selection_1 = data),
                                master$selected_sites_EG))
  }

  # selection process
  if (verbose == TRUE) {
    message("Running algorithm for selecting sites, please wait...")
  }

  while (np > expected_points) {
    # thinning
    thin <- point_thinning(data, variable_1, variable_2, dist, space = "E",
                           max_n_samplings, replicates, set_seed)
    np <- nrow(thin[[1]])
    if (verbose == TRUE) {
      message("    Distance  ", dist, "  resulted in  ", np, "  points")
    }

    if (np <= expected_points) {
      if (np == expected_points) {
        # selection done
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

  # Returning results
  ## adding preselected to selected
  if (use_preselected_sites == TRUE) {
    thin <- lapply(thin, function(x) {
      rbind(master$preselected_sites[, -1], x[, colnames(x) != "Selected_blocks"])
    })
  }

  ## filtering by distances success
  if (length(thin) > 1 & !is.null(median_distance_filter)) {
    thin <- distance_filter(thin, median_distance_filter)
  }

  ## naming and returning
  names(thin) <- paste0("selection_", 1:length(thin))

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
