#' Subset PAM according to selected sites
#'
#' @description
#'
#' @param base_PAM object of class base_PAM obtained using the function
#' \code{\link{base_PAM}}.
#' @param master_selection object of class master_selection. This object can be
#' obtained using the functions: \code{\link{random_selection}},
#' \code{\link{uniformG_selection}}, \code{\link{uniformE_selection}},
#' \code{\link{EG_selection}},
#' @param selection_type type of selection to be considered to subset
#' \code{base_PAM}. Options are: "all", "selected_sites_random", "selected_sites_E",
#' "selected_sites_G", "selected_sites_EG". The default, "all", uses all selection
#' types present in \code{master_selection}.
#'
#' @return
#' An object of class PAM_subset containing the original \code{base_PAM} and
#' other subsets of the PAM according to \code{selection_type}.
#'
#' @usage
#' subset_PAM(base_PAM, master_selection, selection_type = "all")
#'
#' @export
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#' data("species_data", package = "biosurvey")
#'
#' # Create base_pam
#' b_pam <- base_pam(data = species_data, master_matrix = m_matrix, cell_size = 100)
#'
#' # Making blocks for analysis with master matrix
#' m_blocks <- make_blocks(m_matrix, variable_1 = "PC1", variable_2 = "PC2",
#'                         n_cols = 20, n_rows = 20, block_type = "equal_area")
#'
#' # Selections
#' ## Randomly
#' m_selection <- random_selection(m_blocks, n_sites = 20, n_samplings = 2)
#'
#' ## Uniformly in G space
#' m_selection <- uniformG_selection(m_selection, expected_points = 20,
#'                                   max_n_samplings = 1, initial_distance = 165,
#'                                   increase = 1, replicates = 5)
#'
#' ## Uniformly in E space
#' m_selection <- uniformE_selection(m_selection, variable_1 = "PC1", variable_2 = "PC2",
#'                                   selection_from = "block_centroids",
#'                                   expected_points = 15, max_n_samplings = 1,
#'                                   initial_distance = 1, increase = 0.1,
#'                                   replicates = 5)

subset_PAM <- function(base_PAM, master_selection, selection_type = "all") {
  if (missing(base_PAM)) {
    stop("Argument 'base_PAM' must be defined.")
  }
  if (missing(master_selection)) {
    stop("Argument 'master_selection' must be defined.")
  }
  if (!selection_type %in% c("all", "selected_sites_random", "selected_sites_E",
                             "selected_sites_G", "selected_sites_EG")) {
    stop("Argument 'selection_type' is not valid, options are:\n'selected_sites_random', 'selected_sites_E', 'selected_sites_G', or 'selected_sites_EG'.")
  }

  WGS84 <- sp::CRS("+init=epsg:4326")

  if (selection_type == "all") {
    al <- grep()
    for (i in 1:length(al)) {

    }
    bbs_sitepam <- merge(bbs_sites$selected_sites, bbs_pam, by = "routeid")
  }

  if (selection_type == "selected_sites_random") {
    rsel <- master_selection$selected_sites_random
    base_PAM$PAM_selected_sites_random <- lapply(rsel, function(x) {
      xp <- sp::SpatialPointsDataFrame(x[, 1:2], x, proj4string = WGS84)
      xid <- data.frame(ID = base_PAM$PAM[xp, ]@data$ID, x)
      merge(xid, base_PAM$PAM@data, by = "ID")
    })
    names(base_PAM$PAM_selected_sites_random) <- names(rsel)
  }

  if (selection_type == "selected_sites_E") {
    rsel <- master_selection$selected_sites_E
    base_PAM$PAM_selected_sites_E <- lapply(rsel, function(x) {
      xp <- sp::SpatialPointsDataFrame(x[, 1:2], x, proj4string = WGS84)
      xid <- data.frame(ID = base_PAM$PAM[xp, ]@data$ID, x)
      merge(xid, base_PAM$PAM@data, by = "ID")
    })
    names(base_PAM$PAM_selected_sites_E) <- names(rsel)
  }

  if (selection_type == "selected_sites_G") {
    rsel <- master_selection$selected_sites_G
    base_PAM$PAM_selected_sites_G <- lapply(rsel, function(x) {
      xp <- sp::SpatialPointsDataFrame(x[, 1:2], x, proj4string = WGS84)
      xid <- data.frame(ID = base_PAM$PAM[xp, ]@data$ID, x)
      merge(xid, base_PAM$PAM@data, by = "ID")
    })
    names(base_PAM$PAM_selected_sites_G) <- names(rsel)
  }

  if (selection_type == "selected_sites_EG") {
    rsel <- master_selection$selected_sites_EG
    base_PAM$PAM_selected_sites_EG <- lapply(rsel, function(x) {
      xp <- sp::SpatialPointsDataFrame(x[, 1:2], x, proj4string = WGS84)
      xid <- data.frame(ID = base_PAM$PAM[xp, ]@data$ID, x)
      merge(xid, base_PAM$PAM@data, by = "ID")
    })
    names(base_PAM$PAM_selected_sites_EG) <- names(rsel)
  }

  return(structure(base_PAM, class = "PAM_subset"))
}
