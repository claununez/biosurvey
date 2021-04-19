#' Dissimilarity indices from PAM_subset
#'
#' @description computes dissimilarity indices for each set of selected
#' sites contained in elements of \code{PAM_subset} that contain information of
#' species incidence (presence-absence). Calculations are done also among sets
#' of selected sites.
#'
#' @param PAM_subset object of class PAM_subset obtained using the function
#' \code{\link{subset_PAM}}.
#' @param selection_type type of selection to be considered when creating
#' dissimilarity matrices for elements in \code{PAM_subset}. Options are:
#' "all", "random", "E", "G", and "EG". The default, "all", uses all selection
#' types present in \code{PAM_subset}.
#' @param method (character) dissimilarity index to be passed to function
#' \code{\link[vegan]{vegdist}}. Default = "jaccard". See details.
#' @param ... other arguments to be passed to function
#' \code{\link[vegan]{vegdist}}.
#'
#' @return
#' A list containing:
#'
#' - Dissimilarity matrices for all PAMs reduced based on distinct sets of
#' selected sites.
#' - A matrix summarizing incidences from all sets of selected sites.
#' - A dissimilarity matrix for the summary of incidences for all sets of
#' selected sites.
#' - The result of clustering sets of selected sites based on dissimilarities.
#'
#' @details
#' Important details about the process performed to compute dissimilarity
#' indices can be seen in the documentation of \code{\link[vegan]{vegdist}}.
#'
#' @usage
#' selected_sites_DI(PAM_subset, selection_type = "all", method = "jaccard",
#'                   ...)
#'
#' @export
#' @importFrom vegan vegdist
#' @importFrom stats hclust
#'
#' @examples
#' # Data
#' data("b_pam", package = "biosurvey")
#' data("m_selection", package = "biosurvey")
#'
#' # Subsetting base PAM according to selections
#' sub_pam_all <- subset_PAM(b_pam, m_selection, selection_type = "all")
#'
#' # Calculating dissimilarities
#' DI_sel <- selected_sites_DI(sub_pam_all)

selected_sites_DI <- function(PAM_subset, selection_type = "all",
                              method = "jaccard", ...) {

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
  diss <- list()
  dise <- list()
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

  # Dissimilarity calculation
  ## Random
  if ("PAM_selected_sites_random" %in% selection_type &
      !is.null(PAM_subset$PAM_selected_sites_random)) {
    rsel <- PAM_subset$PAM_selected_sites_random
    diss$DI_selected_sites_random  <- lapply(rsel, function(x) {
      comat <- x[, icol:fcol]
      rownames(comat) <- paste0("Site_" , 1:nrow(comat))
      vegan::vegdist(comat, method = method, ...)
    })
    sums <- lapply(rsel, function(x) {colSums(x[, icol:fcol]) > 0})
    sums <- do.call(rbind, sums) * 1
    rownames(sums) <- paste0("random_", rownames(sums))
    dise$random <- sums
  }

  ## G
  if ("PAM_selected_sites_G" %in% selection_type &
      !is.null(PAM_subset$PAM_selected_sites_G)) {
    rsel <- PAM_subset$PAM_selected_sites_G
    diss$DI_selected_sites_G  <- lapply(rsel, function(x) {
      comat <- x[, icol:fcol]
      rownames(comat) <- paste0("Site_" , 1:nrow(comat))
      vegan::vegdist(comat, method = method, ...)
    })
    sums <- lapply(rsel, function(x) {colSums(x[, icol:fcol]) > 0})
    sums <- do.call(rbind, sums) * 1
    rownames(sums) <- paste0("G_", rownames(sums))
    dise$G <- sums
  }

  ## E
  if ("PAM_selected_sites_E" %in% selection_type &
      !is.null(PAM_subset$PAM_selected_sites_E)) {
    rsel <- PAM_subset$PAM_selected_sites_E
    diss$DI_selected_sites_E  <- lapply(rsel, function(x) {
      comat <- x[, icol:fcol]
      rownames(comat) <- paste0("Site_" , 1:nrow(comat))
      vegan::vegdist(comat, method = method, ...)
    })
    sums <- lapply(rsel, function(x) {colSums(x[, icol:fcol]) > 0})
    sums <- do.call(rbind, sums) * 1
    rownames(sums) <- paste0("E_", rownames(sums))
    dise$E <- sums
  }

  ## EG
  if ("PAM_selected_sites_EG" %in% selection_type &
      !is.null(PAM_subset$PAM_selected_sites_EG)) {
    rsel <- PAM_subset$PAM_selected_sites_EG
    diss$DI_selected_sites_EG  <- lapply(rsel, function(x) {
      comat <- x[, icol:fcol]
      rownames(comat) <- paste0("Site_" , 1:nrow(comat))
      vegan::vegdist(comat, method = method, ...)
    })
    sums <- lapply(rsel, function(x) {colSums(x[, icol:fcol]) > 0})
    sums <- do.call(rbind, sums) * 1
    rownames(sums) <- paste0("EG_", rownames(sums))
    dise$EG <- sums
  }

  # dissimilarity among selections
  diss$all_selections <- do.call(rbind, dise)

  ## computing dissimilarities
  diss$DI_selections <- vegan::vegdist(diss$all_selections, method = method,
                                       ...)

  ## clustering according to dissimilarities
  diss$cluster_selections <- stats::hclust(diss$DI_selections)

  return(diss)
}


