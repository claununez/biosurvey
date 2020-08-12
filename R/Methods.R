#' Print a short version of elements in master objects
#' @name print
#' @aliases print,master_matrix-method print,master_selection-method
#' @aliases print,base_PAM-method print,PAM_subset-method
#' @param object object of class master_matrix, master_selection, base_PAM, or
#' PAM_subset.
#' @export
#' @rdname print

print.master_matrix <- function(object) {
  cat("data_matrix:\n")
  print(head(object$data_matrix))
  cat("...\n")

  cat("\npreselected_sites:\n")
  if (!is.null(object$preselected_sites)) {
    print(object$preselected_sites)
  } else {
    cat("Empty\n")
  }

  cat("\nregion:\n")
  print(object$region)

  cat("\nmask:\n")
  if (!is.null(object$mask)) {
    print(object$mask)
  } else {
    cat("Empty\n")
  }

  cat("\nraster_base:\n")
  print(object$raster_base)

  cat("\nPCA_results:\n")
  if (!is.null(object$PCA_results)) {
    print(object$PCA_results)
  } else {
    cat("Empty\n")
  }
}


#' @export
#' @rdname print

print.master_selection <- function(object) {

  print(structure(object[1:6], class = "master_matrix"))

  cat("\nselected_sites_random:\n")
  if (!is.null(object$selected_sites_random)) {
    cat("First of", length(object$selected_sites_random), "element(s).\n")
    print(head(object$selected_sites_random[[1]]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }

  cat("\nselected_sites_G:\n")
  if (!is.null(object$selected_sites_G)) {
    cat("First of", length(object$selected_sites_G), "element(s).\n")
    print(head(object$selected_sites_G[[1]]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }

  cat("\nselected_sites_E:\n")
  if (!is.null(object$selected_sites_E)) {
    cat("First of", length(object$selected_sites_E), "element(s).\n")
    print(head(object$selected_sites_E[[1]]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }

  cat("\nselected_sites_EG:\n")
  if (!is.null(object$selected_sites_EG)) {
    cat("First of", length(object$selected_sites_EG), "element(s).\n")
    print(head(object$selected_sites_EG[[1]]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }
}



#' @export
#' @rdname print

print.base_PAM <- function(object) {
  cat("PAM:\n")
  print(object$PAM)

  cat("\n\nPAM_indices:\n")
  if (!is.null(object$PAM_indices)) {
    cat("  One_value_indices:\n")
    print(na.omit(object$PAM_indices$One_value_indices))

    cat("\n  Richness:\n")
    cat(" ", head(object$PAM_indices$Richness), "...\n")

    cat("\n  Range:\n")
    cat(" ", head(object$PAM_indices$Range), "...\n")

    cat("\n  Richness_normalized:\n")
    cat(" ", head(object$PAM_indices$Richness_normalized), "...\n")

    cat("\n  Range_normalized:\n")
    cat(" ", head(object$PAM_indices$Range_normalized), "...\n")

    cat("\n  Dispersion_field:\n")
    if (!is.null(object$PAM_indices$Dispersion_field)) {
      cat(" ", head(object$PAM_indices$Dispersion_field), "...\n")
    } else {
      cat("  Empty\n")
    }

    cat("\n  Shared_community_composition:\n")
    if (!is.null(object$PAM_indices$Shared_community_composition)) {
      cat(" ", head(object$PAM_indices$Shared_community_composition), "...\n")
    } else {
      cat("  Empty\n")
    }

    cat("\n  Mean_composition_covariance:\n")
    if (!is.null(object$PAM_indices$Mean_composition_covariance)) {
      cat(" ", head(object$PAM_indices$Mean_composition_covariance), "...\n")
    } else {
      cat("  Empty\n")
    }

    cat("\n  Mean_range_covariance:\n")
    if (!is.null(object$PAM_indices$Mean_range_covariance)) {
      cat(" ", head(object$PAM_indices$Mean_range_covariance), "...\n")
    } else {
      cat("  Empty\n")
    }

    cat("\n  Cov_mat_sites_composition:\n")
    if (!is.null(object$PAM_indices$Cov_mat_sites_composition)) {
      print(head(object$PAM_indices$Cov_mat_sites_composition))
      cat("...\n")
    } else {
      cat("  Empty\n")
    }

    cat("\n  Cov_mat_species_ranges\n")
    if (!is.null(object$PAM_indices$Cov_mat_species_ranges)) {
      print(head(object$PAM_indices$Cov_mat_species_ranges))
      cat("...\n")
    } else {
      cat("  Empty\n")
    }

  } else {
    cat("Empty\n")
  }
}


#' @export
#' @rdname print

print.PAM_subset <- function(object) {

  print(structure(object[1:2], class = "base_PAM"))

  cat("\n\nPAM_selected_sites_random:\n")
  if (!is.null(object$PAM_selected_sites_random)) {
    cat("First of", length(object$PAM_selected_sites_random), "element(s).\n")
    print(head(object$PAM_selected_sites_random[[1]][, 1:6]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }

  cat("\nPAM_selected_sites_G:\n")
  if (!is.null(object$PAM_selected_sites_G)) {
    cat("First of", length(object$PAM_selected_sites_G), "element(s).\n")
    print(head(object$PAM_selected_sites_G[[1]][, 1:6]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }

  cat("\nPAM_selected_sites_E:\n")
  if (!is.null(object$PAM_selected_sites_E)) {
    cat("First of", length(object$PAM_selected_sites_E), "element(s).\n")
    print(head(object$PAM_selected_sites_E[[1]][, 1:6]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }

  cat("\nPAM_selected_sites_EG:\n")
  if (!is.null(object$PAM_selected_sites_EG)) {
    cat("First of", length(object$PAM_selected_sites_EG), "element(s).\n")
    print(head(object$PAM_selected_sites_EG[[1]][, 1:6]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }
}


#' Summary of attributes and results
#' @name summary
#' @aliases summary,master_matrix-method summary,master_selection-method
#' @aliases summary,base_PAM-method summary,PAM_subset-method
#' @param object object of class master_matrix or master_selection.
#' @export
#' @return
#' A printed summary.
#' @rdname summary

summary.master_matrix <- function(object) {

  cat("\n                     Summary of a master_matrix object\n")
  cat("---------------------------------------------------------------------------\n\n")
  cat("Data matrix summary:\n")
  print(summary(object$data_matrix))
  if (!is.null(object$preselected_sites)) {
    cat("\n\nSites preselected by user:\n")
    print(object$preselected_sites[, 1:3])
  } else {
    cat("\n\nNo preselected sites were defined\n")
  }
  cat("\n\nRegion of interest:\n")
  print(object$region)
  if (!is.null(object$mask)) {
    cat("\n\nMask to area of interest:\n")
    print(object$mask)
  }
}

#' @export
#' @rdname summary

summary.master_selection <- function(object, nrow = 6, ncol = 2) {

  cat("\n                  Summary of a master_selection object\n")
  cat("---------------------------------------------------------------------------\n\n")
  if (!is.null(object$preselected_sites)) {
    cat("Sites preselected by user:\n")
    print(object$preselected_sites[, 1:3])
  } else {
    cat("No preselected sites were defined\n")
  }

  cat("\n\n")
  cat(ncol, "columns and", nrow, "rows of first elements are shown")

  cat("\n\nSites selected randomly:\n")
  if (!is.null(object$selected_sites_random)) {
    print(object$selected_sites_random[[1]][1:nrow, 1:ncol])
  } else {
    cat("Empty\n")
  }

  cat("\n\nSites selected uniformly in G space:\n")
  if (!is.null(object$selected_sites_G)) {
    print(object$selected_sites_G[[1]][1:nrow, 1:ncol])
  } else {
    cat("Empty\n")
  }

  cat("\n\nSites selected uniformly in E space:\n")
  if (!is.null(object$selected_sites_E)) {
    print(object$selected_sites_E[[1]][1:nrow, 1:ncol])
  } else {
    cat("Empty\n")
  }

  cat("\n\nSites selected uniformly in E space, considering G structure:\n")
  if (!is.null(object$selected_sites_EG)) {
    print(object$selected_sites_EG[[1]][1:nrow, 1:ncol])
  } else {
    cat("Empty\n")
  }
}



#' @export
#' @rdname summary

summary.base_PAM <- function(object) {

  cat("\n                      Summary of a base_PAM object\n")
  cat("---------------------------------------------------------------------------\n\n")
  cat("Presence-absence matrix:\n")
  cat("  Number of cells:   ", object$PAM_indices$One_value_indices["Sites_Cells", 1])
  cat("\n  Number of species: ", object$PAM_indices$One_value_indices["Species", 1])

  cat("\n\nSpatial object representing the PAM:\n")
  print(object$PAM)
}


#' @export
#' @rdname summary

summary.PAM_subset <- function(object) {

  cat("\n                      Summary of a PAM_subset object\n")
  cat("---------------------------------------------------------------------------\n\n")
  cat("Complete presence-absence matrix:\n")
  cat("  Number of cells:   ", object$PAM_indices$One_value_indices["Sites_Cells", 1])
  cat("\n  Number of species: ", object$PAM_indices$One_value_indices["Species", 1])

  cat("\n\nSubsets of PAM:\n")
  sele <- object[3:6]
  snames <- names(sele)
  nnull <- which(!sapply(sele, is.null))
  ncells <- sapply(sele[nnull], function(x) {nrow(x[[1]])})

  sps <- sapply(sele[nnull], function(x) {
    cnam <- colnames(x[[1]])
    icol <- which(cnam == "Latitude_PAM") + 1
    cnam <- cnam[icol:length(cnam)]
    sum(apply(x[[1]][, icol:ncol(x[[1]])], 2, max) == 1)
  })
  print(data.frame(Cells = ncells, Species = sps))
}
