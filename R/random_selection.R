#' Random selection of survey sites
#'
#' @description Random selection of sites to be sampled in a survey. Sites are
#' selected from a set of points provided in \code{master}.
#'
#' @param master a master_matrix object derived from function
#' \code{\link{master_matrix}} or a master_selection object derived from functions
#' \code{\link{uniformG_selection}}, \code{\link{uniformE_selection}}
#' or \code{\link{EG_selection}}.
#' @param n_sites (numeric) number of sites to be selected from
#' \code{master_matrix} to be used as sites to be sampled in survey.
#' @param n_samplings (numeric) number of processes of selection, which will
#' turn into multiple options for a process of survey planning. Default = 1.
#' @param median_distance_filter (character) optional argument to define a median
#' distance-based filter based on which sets of sampling sites will be selected.
#' The default, NULL, does not apply such a filter. Options are: "max" and "min".
#' @param set_seed (numeric) integer value to specify a initial seed. Default = 1.
#'
#' @return
#' A master_selection object (S3) with an additional element called
#' selected_sites_random containing one or more sets of selected sites.
#'
#' @details
#' Survey sites are selected randomly from the entire set of points provided in
#' \code{master$master_matrix}. Considering the environmental space, as points
#' are selected randomly, sites that have environmental conditions that are common
#' in the region of interest will be sampled more that other sites that present
#' condition that are not as common.
#'
#' To see how common or rare are distinct environments in the region of interest
#' the function \code{\link{explore_data_EG}} can be used. Common environmental
#' conditions are those that are present in areas of higher density in one of the
#' plots obtained with \code{\link{explore_data_EG}}.
#'
#' @seealso
#' \code{\link{uniformG_selection}}, \code{\link{uniformE_selection}},
#' \code{\link{EG_selection}}, \code{\link{plot_sites_EG}}
#'
#' @usage
#' random_selection(master, n_sites, n_samplings = 1,
#'                  median_distance_filter = NULL, set_seed = 1)
#'
#' @export
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#'
#' r_selection <- random_selection(m_matrix, n_sites = 20, n_samplings = 5)

random_selection <- function(master, n_sites, n_samplings = 1,
                             median_distance_filter = NULL, set_seed = 1) {

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

  # Selection of sites
  data <- master$master_matrix

  selected_sites <- lapply(1:n_samplings, function(x) {
    set.seed(set_seed + x - 1)
    sam <- sample(nrow(data), n_sites)
    dat <- data[sam, ]
  })
  names(selected_sites) <- paste0("selection_", 1:n_samplings)

  # Post filtering of sites according to distance argument
  if (length(selected_sites) > 1 & !is.null(median_distance_filter)) {
    selected_sites <- distance_filter(selected_sites, median_distance_filter)

  }

  # Returning results
  master$selected_sites_random <- selected_sites

  return(structure(master, class = "master_selection"))
}
