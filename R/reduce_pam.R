#' Subset PAM according to selected sites
#'
#' @description Subsets of a base_PAM object according to survey sites contained
#' in a master_selection object.
#'
#' @param base_PAM object of class base_PAM obtained using the function
#' \code{\link{base_PAM}}.
#' @param master_selection object of class master_selection. This object can be
#' obtained using the functions: \code{\link{random_selection}},
#' \code{\link{uniformG_selection}}, \code{\link{uniformE_selection}}, and
#' \code{\link{EG_selection}}.
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
#' data("b_pam", package = "biosurvey")
#' data("m_selection", package = "biosurvey")
#'
#' # Subsetting base PAM according to selections
#' ## only uniform in G
#' sub_pam_G <- subset_PAM(b_pam, m_selection, selection_type = "selected_sites_G")
#'
#' ## all at the time
#' sub_pam_all <- subset_PAM(b_pam, m_selection, selection_type = "all")

subset_PAM <- function(base_PAM, master_selection, selection_type = "all") {

  # Initial tests
  if (missing(base_PAM)) {
    stop("Argument 'base_PAM' must be defined.")
  }
  if (missing(master_selection)) {
    stop("Argument 'master_selection' must be defined.")
  }
  if (class(master_selection)[1] != "master_selection") {
    stop("Object 'master_selection' must be of class 'master_selection'.")
  }
  if (!selection_type %in% c("all", "selected_sites_random", "selected_sites_E",
                             "selected_sites_G", "selected_sites_EG")) {
    stop("Argument 'selection_type' is not valid, options are:\n'selected_sites_random', 'selected_sites_E', 'selected_sites_G', or 'selected_sites_EG'.")
  }

  # Identifying selection types if all
  if (selection_type == "all") {
    selects <- names(master_selection)
    selection_type <- grep("selected_sites", selects, value = TRUE)
  }

  # Joining PAM with selected sites
  ## random
  if ("selected_sites_random" %in% selection_type) {
    rpsel <- master_selection$selected_sites_random
    rpsel <- selected_sites_PAM(rpsel, base_PAM)
  } else {
    rpsel <- NULL
  }

  ## E
  if ("selected_sites_G" %in% selection_type) {
    gpsel <- master_selection$selected_sites_E
    gpsel <- selected_sites_PAM(gpsel, base_PAM)
  } else {
    gpsel <- NULL
  }

  ## G
  if ("selected_sites_E" %in% selection_type) {
    epsel <- master_selection$selected_sites_G
    epsel <- selected_sites_PAM(epsel, base_PAM)
  } else {
    epsel <- NULL
  }

  ## EG
  if ("selected_sites_EG" %in% selection_type) {
    egpsel <- master_selection$selected_sites_EG
    egpsel <- selected_sites_PAM(egpsel, base_PAM)
  } else {
    egpsel <- NULL
  }

  # Returning results
  return(new_PAM_subset(base_PAM$PAM, base_PAM$PAM_indices, rpsel, gpsel,
                        epsel, egpsel))
}
