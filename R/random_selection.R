#' Random selection of survey sites
#'
#' @description Random selection of sites to be sampled in a survey. Sites are
#' selected from a set of points provided in \code{master}.
#'
#' @param master a master_matrix object derived from function
#' \code{\link{prepare_master_matrix}} or a master_selection object derived
#' from functions \code{\link{uniformG_selection}},
#' \code{\link{uniformE_selection}} or \code{\link{EG_selection}}.
#' @param n_sites (numeric) number of sites to be selected from
#' \code{master_matrix} to be used as sites to be sampled in survey.
#' @param n_samplings (numeric) number of processes of selection, which will
#' turn into multiple options for a process of survey planning. Default = 1.
#' @param use_preselected_sites (logical) whether to use sites that have been
#' defined as part of the selected sites previous any selection. Object in
#' \code{master} must contain the site(s) preselected in and element of name
#' "preselected_sites" for this argument to be effective. Default = TRUE.
#' @param median_distance_filter (character) optional argument to define a
#' median distance-based filter based on which sets of sampling sites will be
#' selected. The default, NULL, does not apply such a filter. Options are:
#' "max" and "min".
#' @param set_seed (numeric) integer value to specify a initial seed.
#' Default = 1.
#' @param verbose (logical) whether or not to print messages about the process.
#' Default = TRUE.
#'
#' @return
#' A master_selection object (S3) with an additional element called
#' selected_sites_random containing one or more sets of selected sites.
#'
#' @details
#' Survey sites are selected randomly from the entire set of points provided in
#' \code{master$data_matrix}. Considering the environmental space, as points
#' are selected randomly, sites that have environmental conditions that are
#' common in the region of interest will be sampled more that other sites that
#' present condition that are not as common.
#'
#' To see how common or rare are distinct environments in the region of interest
#' the function \code{\link{explore_data_EG}} can be used. Common environmental
#' conditions are those that are present in areas of higher density in one of
#' the plots obtained with \code{\link{explore_data_EG}}.
#'
#' As multiple sets could result from selection, the argument of the function
#' \code{median_distance_filter} could be used to select the set of sites with
#' the maximum ("max") or minimum ("min") median distance among selected sites.
#' Option "max" will increase the geographic distance among sampling sites,
#' which could be desirable if the goal is to cover the region of interest more
#' broadly. The other option "min", could be used in cases when the goal is to
#' reduce resources and time needed to sample such sites.
#'
#' @seealso
#' \code{\link{uniformG_selection}}, \code{\link{uniformE_selection}},
#' \code{\link{EG_selection}}, \code{\link{plot_sites_EG}}
#'
#' @usage
#' random_selection(master, n_sites, n_samplings = 1,
#'                  use_preselected_sites = TRUE, median_distance_filter = NULL,
#'                  set_seed = 1, verbose = TRUE)
#'
#' @export
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#'
#' r_selection <- random_selection(m_matrix, n_sites = 20, n_samplings = 5)

random_selection <- function(master, n_sites, n_samplings = 1,
                             use_preselected_sites = TRUE,
                             median_distance_filter = NULL,
                             set_seed = 1, verbose = TRUE) {

  # Initial tests
  if (missing(master)) {
    stop("Argument 'master' must be defined.")
  }
  if (missing(n_sites)) {
    stop("Argument 'n_sites' must be defined.")
  }
  if (!is.null(median_distance_filter)) {
    if (!median_distance_filter %in% c("max", "min")) {
      stop("Argument 'median_distance_filter' is not valid, see function's help.")
    }
  }
  if (use_preselected_sites == TRUE & is.null(master$preselected_sites)) {
    message("Element 'preselected_sites' in 'master' is NULL, setting\n'use_preselected_sites' = FALSE")
    use_preselected_sites <- FALSE
  }

  # Selection of sites
  if (verbose == TRUE) {
    message("Selecting sampling sites randomly")
  }

  data <- master$data_matrix
  n <- nrow(data)

  if (use_preselected_sites == TRUE) {
    pre <- master$preselected_sites
  }

  selected_sites <- lapply(1:n_samplings, function(x) {
    set.seed(set_seed + x - 1)
    sam <- sample(n, n_sites)
    if (use_preselected_sites == TRUE) {
      dat <- unique(rbind(pre[, -1], data[sam, ]))
      dat <- dat[1:n_sites, ]
    } else {
      dat <- data[sam, ]
    }
  })

  # Post filtering of sites according to distance argument
  if (length(selected_sites) > 1 & !is.null(median_distance_filter)) {
    selected_sites <- distance_filter(selected_sites, median_distance_filter)

  }

  # Returning results
  names(selected_sites) <- paste0("selection_", 1:length(selected_sites))

  ## arguments as attributes
  other_args <- list(arguments = list(n_sites = n_sites,
                                      n_samplings = n_samplings,
                                      use_preselected_sites = use_preselected_sites,
                                      median_distance_filter = median_distance_filter,
                                      set_seed = set_seed))
  attributes(selected_sites) <- c(attributes(selected_sites), other_args)

  if (verbose == TRUE) {
    message("Total number of sites selected: ", nrow(selected_sites[[1]]))
  }

  return(new_master_selection(master$data_matrix, master$preselected_sites,
                              master$region, master$mask, master$raster_base,
                              master$PCA_results, selected_sites,
                              master$selected_sites_G, master$selected_sites_E,
                              master$selected_sites_EG))
}
