#' Random selection of survey sites
#'
#' @description Random selection of sites to be sampled in a survey. Sites are
#' selected from a set of points that are provided in one of the arguments.
#'
#' @param master_matrix object derived from function \code{\link{master_matrix}}.
#' Optionally, if master_matrix is not necessary, a list containing an object of
#' class data.frame with at least two columns represening two variables. The name
#' of this element in the list must be "master_matrix". For instance:
#' \code{my_list <- list(master_matrix = YOUR_data.frame)}.
#' @param n_sites (numeric) number of sites to be selected from
#' \code{master_matrix} to be used as sites to be sampled in survey.
#' @param n_samplings (numeric) number of processes of selection, which will
#' turn into multiple options for a process of survey planning. Default = 1.
#' @param set_seed (numeric) integer value to specify a initial seed. Default = 1.
#'
#' @return
#' The master_matrix list with an aditional element containing one or more sets
#' of selected sites.
#'
#' @usage
#' random_selection(master_matrix, n_sites, n_samplings = 1, set_seed = 1)
#'
#' @export
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#'
#' r_selection <- random_selection(m_matrix, n_sites = 20, n_samplings = 5)

random_selection <- function(master_matrix, n_sites, n_samplings = 1,
                             set_seed = 1) {
  if (missing(master_matrix)) {
    stop("Argument 'master_matrix' must be defined.")
  }
  if (missing(n_sites)) {
    stop("Argument 'n_sites' must be defined.")
  }
  data <- master_matrix$master_matrix

  selected_sites <- lapply(1:n_samplings, function(x) {
    set.seed(set_seed + x - 1)
    sam <- sample(nrow(data), n_sites)
    dat <- data[sam, ]
  })
  names(selected_sites) <- paste0("selection_", 1:n_samplings)

  master_matrix$selected_sites <- selected_sites

  return(master_matrix)
}
