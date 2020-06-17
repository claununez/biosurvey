#' Random selection of survey sites
#'
#' @description Random selection of sites to be sampled in a survey. Sites are
#' selected from a set of points that are provided in one of the arguments.
#'
#' @param master a master_matrix object derived from function
#' \code{\link{master_matrix}} or a master_selection object derived from functions
#' \code{\link{uniformG_selection}}, \code{\link{uniformE_selection}}
#' or \code{\link{EG_selection}}.
#' @param n_sites (numeric) number of sites to be selected from
#' \code{master_matrix} to be used as sites to be sampled in survey.
#' @param n_samplings (numeric) number of processes of selection, which will
#' turn into multiple options for a process of survey planning. Default = 1.
#' @param set_seed (numeric) integer value to specify a initial seed. Default = 1.
#'
#' @return
#' A master_selection object (S3) with an aditional element called
#' selected_sites_random containing one or more sets of selected sites.
#'
#' @usage
#' random_selection(master, n_sites, n_samplings = 1, set_seed = 1)
#'
#' @export
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#'
#' r_selection <- random_selection(m_matrix, n_sites = 20, n_samplings = 5)

random_selection <- function(master, n_sites, n_samplings = 1, set_seed = 1) {

  # Initial tests
  if (missing(master)) {
    stop("Argument 'master' must be defined.")
  }
  if (missing(n_sites)) {
    stop("Argument 'n_sites' must be defined.")
  }
  data <- master$master_matrix

  selected_sites <- lapply(1:n_samplings, function(x) {
    set.seed(set_seed + x - 1)
    sam <- sample(nrow(data), n_sites)
    dat <- data[sam, ]
  })
  names(selected_sites) <- paste0("selection_", 1:n_samplings)

  master$selected_sites_random <- selected_sites

  return(structure(master, class = "master_selection"))
}
