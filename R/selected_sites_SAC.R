#' Species accumulation curves from PAM_subset
#'
#' @description Creates species accumulation curves for each set of selected sites
#' contained in elements of \code{PAM_subset} that contain information of
#' species incidence (presence-absence).
#'
#' @param PAM_subset object of class PAM_subset obtained using the function
#' \code{\link{subset_PAM}}.
#' @param selection_type type of selection to be considered when creating SAC for
#' elements in \code{PAM_subset}. Options are: "all", "random", "E", "G", "EG".
#' The default, "all", uses all selection types present in
#' \code{PAM_subset}.
#' @param method (character) species accumulation method to be passed to function
#' \code{\link[vegan]{specaccum}}. Default = "exact".
#' @param ... other arguments to be passed to function \code{\link[vegan]{specaccum}}.
#'
#' @return
#' A list of species accumulation curves (SACs, "\code{specaccum}" objects) for
#' all sets of selected sites according to option defined in \code{selection_type}.
#'
#' @details
#' Important details about the process performed to obtain each of the SACs can
#' be seen in the help for function \code{\link[vegan]{specaccum}}.
#'
#' @usage
#' selected_sites_SAC(PAM_subset, selection_type = "all", method = "exact", ...)
#'
#' @export
#' @importFrom vegan specaccum
#'
#' @examples
#' # Data
#' data("b_pam", package = "biosurvey")
#' data("m_selection", package = "biosurvey")
#'
#' # Subsetting base PAM according to selections
#' sub_pam_all <- subset_PAM(b_pam, m_selection, selection_type = "all")
#'
#' SACs <- selected_sites_SAC(PAM_subset = sub_pam_all, selection_type = "all")

selected_sites_SAC <- function(PAM_subset, selection_type = "all",
                               method = "exact", ...) {
  # Initial tests
  if (missing(PAM_subset)) {
    stop("Argument 'PAM_subset' must be defined.")
  }
  if (class(PAM_subset)[1] != "PAM_subset") {
    stop("Object 'PAM_subset' must be of class 'PAM_subset'.")
  }
  if (!selection_type %in% c("all", "random", "E", "G", "EG")) {
    stop("Argument 'selection_type' is not valid, options are: 'all'', 'random', 'E', 'G', or 'EG'.")
  } else {
    if (!selection_type == "all") {
      selection_type <- paste0("selected_sites_", selection_type)
    }
  }

  # Initial pre-processing
  sac <- list()
  nnull <- which(!sapply(PAM_subset[3:6], is.null))[1] + 2

  rsel <- PAM_subset[[nnull]]
  fcol <- ncol(rsel[[1]])
  icol <- which(colnames(rsel[[1]]) == "Latitude_PAM") + 1

  # Identifying selection types if all
  if (selection_type == "all") {
    selects <- names(PAM_subset)
    selection_type <- grep("selected_sites", selects, value = TRUE)
  } else {
    selection_type <- paste0("PAM_", selection_type)
  }

  # SAC calculations
  ## random
  if ("PAM_selected_sites_random" %in% selection_type &
      !is.null(PAM_subset$PAM_selected_sites_random)) {
    rsel <- PAM_subset$PAM_selected_sites_random
    sac$SAC_selected_sites_random  <- lapply(rsel, function(x) {
      vegan::specaccum(comm = x[, icol:fcol], method = method, ...)
    })
  }

  ## E
  if ("PAM_selected_sites_E" %in% selection_type &
      !is.null(PAM_subset$PAM_selected_sites_E)) {
    rsel <- PAM_subset$PAM_selected_sites_E
    sac$SAC_selected_sites_E  <- lapply(rsel, function(x) {
      vegan::specaccum(comm = x[, icol:fcol], method = method, ...)
    })
  }

  ## G
  if ("PAM_selected_sites_G" %in% selection_type &
      !is.null(PAM_subset$PAM_selected_sites_G)) {
    rsel <- PAM_subset$PAM_selected_sites_G
    sac$SAC_selected_sites_G  <- lapply(rsel, function(x) {
      vegan::specaccum(comm = x[, icol:fcol], method = method, ...)
    })
  }

  ## EG
  if ("PAM_selected_sites_EG" %in% selection_type &
      !is.null(PAM_subset$PAM_selected_sites_EG)) {
    rsel <- PAM_subset$PAM_selected_sites_EG
    sac$SAC_selected_sites_EG  <- lapply(rsel, function(x) {
      vegan::specaccum(comm = x[, icol:fcol], method = method, ...)
    })
  }

  return(sac)
}
