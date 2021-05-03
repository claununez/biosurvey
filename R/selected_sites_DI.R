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
#' @param verbose (logical) whether or not to print messages about the process.
#' Default = TRUE.
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
#'                   verbose = TRUE, ...)
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
                              method = "jaccard", verbose = TRUE, ...) {

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
  if (verbose == TRUE) {
    message("Running analysis...")
  }
  ## Random
  if ("PAM_selected_sites_random" %in% selection_type &
      !is.null(PAM_subset$PAM_selected_sites_random)) {
    if (verbose == TRUE) {
      message("Random selection")
    }
    rsel <- PAM_subset$PAM_selected_sites_random
    diss$DI_selected_sites_random <- dis_loop(rsel, icol, fcol, method, verbose,
                                              ...)
    diss$cluster_random <- lapply(diss$DI_selected_sites_random, function(x) {
      stats::hclust(x)
    })
    sums <- lapply(rsel, function(x) {colSums(x[, icol:fcol]) > 0})
    sums <- do.call(rbind, sums) * 1
    rownames(sums) <- paste0("random_", rownames(sums))
    dise$random <- sums
  }

  ## G
  if ("PAM_selected_sites_G" %in% selection_type &
      !is.null(PAM_subset$PAM_selected_sites_G)) {
    if (verbose == TRUE) {
      message("G selection")
    }
    rsel <- PAM_subset$PAM_selected_sites_G
    diss$DI_selected_sites_G  <- dis_loop(rsel, icol, fcol, method, verbose,
                                          ...)
    diss$cluster_G <- lapply(diss$DI_selected_sites_G, function(x) {
      stats::hclust(x)
    })
    sums <- lapply(rsel, function(x) {colSums(x[, icol:fcol]) > 0})
    sums <- do.call(rbind, sums) * 1
    rownames(sums) <- paste0("G_", rownames(sums))
    dise$G <- sums
  }

  ## E
  if ("PAM_selected_sites_E" %in% selection_type &
      !is.null(PAM_subset$PAM_selected_sites_E)) {
    if (verbose == TRUE) {
      message("E selection")
    }
    rsel <- PAM_subset$PAM_selected_sites_E
    diss$DI_selected_sites_E  <- dis_loop(rsel, icol, fcol, method, verbose,
                                          ...)
    diss$cluster_E <- lapply(diss$DI_selected_sites_E, function(x) {
      stats::hclust(x)
    })
    sums <- lapply(rsel, function(x) {colSums(x[, icol:fcol]) > 0})
    sums <- do.call(rbind, sums) * 1
    rownames(sums) <- paste0("E_", rownames(sums))
    dise$E <- sums
  }

  ## EG
  if ("PAM_selected_sites_EG" %in% selection_type &
      !is.null(PAM_subset$PAM_selected_sites_EG)) {
    if (verbose == TRUE) {
      message("EG selection")
    }
    rsel <- PAM_subset$PAM_selected_sites_EG
    diss$DI_selected_sites_EG  <- dis_loop(rsel, icol, fcol, method, verbose,
                                           ...)
    diss$cluster_EG <- lapply(diss$DI_selected_sites_EG, function(x) {
      stats::hclust(x)
    })
    sums <- lapply(rsel, function(x) {colSums(x[, icol:fcol]) > 0})
    sums <- do.call(rbind, sums) * 1
    rownames(sums) <- paste0("EG_", rownames(sums))
    dise$EG <- sums
  }

  # dissimilarity among selections
  if (verbose == TRUE) {
    message("Summary of all selections")
  }

  diss$all_selections <- do.call(rbind, dise)

  toex <- which(apply(diss$all_selections, 1, sum) == 0)
  if (length(toex) > 0) {
    if (verbose == TRUE) {
      message("\tOne or more sites were excluded due to lack of species data:\n\t",
              paste(names(toex), collapse = ", "))
    }
    diss$all_selections <- diss$all_selections[-toex, ]
  }

  ## computing dissimilarities
  diss$DI_selections <- vegan::vegdist(diss$all_selections, method = method,
                                       ...)

  ## clustering according to dissimilarities
  diss$cluster_selections <- stats::hclust(diss$DI_selections)

  return(diss)
}


#' Helper to calculate dissimilarities in loop
#' @param site_spp_list list of presence absence matrices for a set of sites or
#' cells of a grid.
#' @param icol number of column where the species list (other columns) starts.
#' @param fcol number of column where the species list (other columns) ends.
#' @param method (character) dissimilarity index to be passed to function
#' \code{\link[vegan]{vegdist}}. Default = "jaccard". See details.
#' @param verbose (logical) whether or not to print messages about the process.
#' Default = TRUE.
#' @param ... other arguments to be passed to function
#' \code{\link[vegan]{vegdist}}.
#'
#' @return A list of results from \code{\link[vegan]{vegdist}}.
#' @usage
#' dis_loop(site_spp_list, icol, fcol, method = "jaccard", verbose = TRUE, ...)
#'
#' @export
#' @importFrom vegan vegdist

dis_loop <- function(site_spp_list, icol, fcol, method = "jaccard",
                     verbose = TRUE, ...) {
  if (missing(site_spp_list)) {stop("Argument 'site_spp_list' is missing")}

  lapply(site_spp_list, function(x) {
    comat <- x[, icol:fcol]
    rownames(comat) <- paste0("Site_" , 1:nrow(comat))
    toex <- which(apply(comat, 1, sum) == 0)
    if (length(toex) > 0) {
      if (verbose == TRUE) {
        message("\tOne or more sites were excluded due to lack of species data:\n\t",
                paste(names(toex), collapse = ", "))
      }
      comat <- comat[-toex, ]
    }
    vegan::vegdist(comat, method = method, ...)
  })
}
