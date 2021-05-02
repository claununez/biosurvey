#' Print a short version of elements in master objects
#' @name print
#' @aliases print,master_matrix-method print,master_selection-method
#' @aliases print,base_PAM-method print,PAM_subset-method
#' @param x object of class master_matrix, master_selection, base_PAM, or
#' PAM_subset.
#' @param ... further arguments to be passed to or from other methods. Ignored
#' in these functions.
#' @export
#' @importFrom methods new
#' @importFrom utils head
#' @rdname print

print.master_matrix <- function(x, ...) {
  if (missing(x)) {"Argument 'x' is missing"}

  cat("data_matrix:\n")
  print(head(x$data_matrix))
  cat("...\n")

  cat("\npreselected_sites:\n")
  if (!is.null(x$preselected_sites)) {
    print(x$preselected_sites)
  } else {
    cat("Empty\n")
  }

  cat("\nregion:\n")
  print(x$region)

  cat("\nmask:\n")
  if (!is.null(x$mask)) {
    print(x$mask)
  } else {
    cat("Empty\n")
  }

  cat("\nraster_base:\n")
  print(x$raster_base)

  cat("\nPCA_results:\n")
  if (!is.null(x$PCA_results)) {
    print(x$PCA_results)
  } else {
    cat("Empty\n")
  }
}


#' @export
#' @rdname print

print.master_selection <- function(x, ...) {
  if (missing(x)) {"Argument 'x' is missing"}

  nc <- ncol(x)
  nc <- ifelse(nc > 6, 6, nc)
  print(structure(x[1:nc], class = "master_matrix"))

  cat("\nselected_sites_random:\n")
  if (!is.null(x$selected_sites_random)) {
    cat("First of", length(x$selected_sites_random), "element(s).\n")
    print(head(x$selected_sites_random[[1]]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }

  cat("\nselected_sites_G:\n")
  if (!is.null(x$selected_sites_G)) {
    cat("First of", length(x$selected_sites_G), "element(s).\n")
    print(head(x$selected_sites_G[[1]]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }

  cat("\nselected_sites_E:\n")
  if (!is.null(x$selected_sites_E)) {
    cat("First of", length(x$selected_sites_E), "element(s).\n")
    print(head(x$selected_sites_E[[1]]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }

  cat("\nselected_sites_EG:\n")
  if (!is.null(x$selected_sites_EG)) {
    cat("First of", length(x$selected_sites_EG), "element(s).\n")
    print(head(x$selected_sites_EG[[1]]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }
}



#' @export
#' @rdname print

print.base_PAM <- function(x, ...) {
  if (missing(x)) {"Argument 'x' is missing"}

  cat("PAM:\n")
  print(x$PAM)

  cat("\n\nPAM_indices:\n")
  if (!is.null(x$PAM_indices)) {
    cat("  One_value_indices:\n")
    print(na.omit(x$PAM_indices$One_value_indices))

    cat("\n  Richness:\n")
    cat(" ", head(x$PAM_indices$Richness), "...\n")

    cat("\n  Range:\n")
    cat(" ", head(x$PAM_indices$Range), "...\n")

    cat("\n  Richness_normalized:\n")
    cat(" ", head(x$PAM_indices$Richness_normalized), "...\n")

    cat("\n  Range_normalized:\n")
    cat(" ", head(x$PAM_indices$Range_normalized), "...\n")

    cat("\n  Dispersion_field:\n")
    if (!is.null(x$PAM_indices$Dispersion_field)) {
      cat(" ", head(x$PAM_indices$Dispersion_field), "...\n")
    } else {
      cat("  Empty\n")
    }

    cat("\n  Shared_community_composition:\n")
    if (!is.null(x$PAM_indices$Shared_community_composition)) {
      cat(" ", head(x$PAM_indices$Shared_community_composition), "...\n")
    } else {
      cat("  Empty\n")
    }

    cat("\n  Mean_composition_covariance:\n")
    if (!is.null(x$PAM_indices$Mean_composition_covariance)) {
      cat(" ", head(x$PAM_indices$Mean_composition_covariance), "...\n")
    } else {
      cat("  Empty\n")
    }

    cat("\n  Mean_range_covariance:\n")
    if (!is.null(x$PAM_indices$Mean_range_covariance)) {
      cat(" ", head(x$PAM_indices$Mean_range_covariance), "...\n")
    } else {
      cat("  Empty\n")
    }

    cat("\n  Cov_mat_sites_composition:\n")
    if (!is.null(x$PAM_indices$Cov_mat_sites_composition)) {
      nc <- ncol(x$PAM_indices$Cov_mat_sites_composition)
      nc <- ifelse(nc > 6, 6, nc)
      print(head(x$PAM_indices$Cov_mat_sites_composition[, 1:nc]))
      cat("...\n")
    } else {
      cat("  Empty\n")
    }

    cat("\n  Cov_mat_species_ranges\n")
    if (!is.null(x$PAM_indices$Cov_mat_species_ranges)) {
      nc <- ncol(x$PAM_indices$Cov_mat_species_ranges)
      nc <- ifelse(nc > 6, 6, nc)
      print(head(x$PAM_indices$Cov_mat_species_ranges[, 1:nc]))
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

print.PAM_subset <- function(x, ...) {
  if (missing(x)) {"Argument 'x' is missing"}

  print(structure(x[1:2], class = "base_PAM"))

  cat("\n\nPAM_selected_sites_random:\n")
  if (!is.null(x$PAM_selected_sites_random)) {
    cat("First of", length(x$PAM_selected_sites_random), "element(s).\n")
    nc <- ncol(x$PAM_selected_sites_random[[1]])
    nc <- ifelse(nc > 6, 6, nc)
    print(head(x$PAM_selected_sites_random[[1]][, 1:nc]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }

  cat("\nPAM_selected_sites_G:\n")
  if (!is.null(x$PAM_selected_sites_G)) {
    cat("First of", length(x$PAM_selected_sites_G), "element(s).\n")
    nc <- ncol(x$PAM_selected_sites_G[[1]])
    nc <- ifelse(nc > 6, 6, nc)
    print(head(x$PAM_selected_sites_G[[1]][, 1:nc]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }

  cat("\nPAM_selected_sites_E:\n")
  if (!is.null(x$PAM_selected_sites_E)) {
    cat("First of", length(x$PAM_selected_sites_E), "element(s).\n")
    nc <- ncol(x$PAM_selected_sites_E[[1]])
    nc <- ifelse(nc > 6, 6, nc)
    print(head(x$PAM_selected_sites_E[[1]][, 1:nc]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }

  cat("\nPAM_selected_sites_EG:\n")
  if (!is.null(x$PAM_selected_sites_EG)) {
    cat("First of", length(x$PAM_selected_sites_EG), "element(s).\n")
    nc <- ncol(x$PAM_selected_sites_EG[[1]])
    nc <- ifelse(nc > 6, 6, nc)
    print(head(x$PAM_selected_sites_EG[[1]][, 1:nc]))
    cat("...\n")
  } else {
    cat("Empty\n")
  }
}


#' @export
#' @rdname print

print.PAM_CS <- function(x, ...) {
  if (missing(x)) {"Argument 'x' is missing"}

  cat("Species:  ", x$Species)
  cat("\nSites_cells:  ", x$Sites_cells)
  cat("\nBeta_W:  ", x$Sites_cells)
  cat("\nSpearman_cor:  ", x$Sites_cells)

  cat("\n\nTheoretical_boundaries:\n")
  cat(" X:\n")
  cat(" ", x$Theoretical_boundaries$x)
  cat("\n Y:\n")
  cat(" ", x$Theoretical_boundaries$y)

  cat("\n\nRichness_normalized:\n")
  cat(head(x$Richness_normalized), "...\n")

  cat("\nDispersion_field_normalized:\n")
  cat(head(x$Dispersion_field_normalized), "...\n")

  cat("\nS_significance_id:\n")
  if (!all(is.na(x$S_significance_id))) {
    cat(head(x$S_significance_id), "...\n")
  } else {
    cat("Empty\n")
  }

  cat("\nRandomized_DF:\n")
  if (!all(is.na(c(x$Randomized_DF)))) {
    nc <- ncol(x$Randomized_DF)
    nc <- ifelse(nc > 6, 6, nc)
    print(head(x$Randomized_DF[, 1:nc]))
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
#' @param nrow number of rows to be printed for selected_sites in a
#' master_selection object.
#' @param ncol number of columns to be printed for selected_sites in a
#' master_selection object.
#' @param ... additional arguments affecting the summary produced. Ignored in
#' these functions.
#' @export
#' @return
#' A printed summary.
#' @rdname summary

summary.master_matrix <- function(object, ...) {
  if (missing(object)) {"Argument 'object' is missing"}

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

summary.master_selection <- function(object, nrow = 6, ncol = 2, ...) {
  if (missing(object)) {"Argument 'object' is missing"}

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

summary.base_PAM <- function(object, ...) {
  if (missing(object)) {"Argument 'object' is missing"}

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

summary.PAM_subset <- function(object, ...) {
  if (missing(object)) {"Argument 'object' is missing"}

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
  print(data.frame(Cells = ncells, Species = sps), row.names = FALSE)
}



#' @export
#' @rdname summary

summary.PAM_CS <- function(object, ...) {
  if (missing(object)) {"Argument 'object' is missing"}

  cat("\n                     Summary of a PAM_CS object\n")
  cat("---------------------------------------------------------------------------\n\n")
  cat("Descriptive values:\n")
  cat("  Number of species:  ", object$Species)
  cat("\n  Number of cells:  ", object$Sites_cells)
  cat("\n  Whittaker's beta:  ", object$Sites_cells)
  cat("\n  Spearman's correlation:  ", object$Sites_cells)

  cat("\n\nBoundaries:\n")
  print(data.frame(x = object$Theoretical_boundaries$x,
                   y = object$Theoretical_boundaries$y), row.names = FALSE)

  cat("\nSummary normalized richness:\n")
  print(summary(object$Richness_normalized))

  cat("\nSummary normalized dispersion field:\n")
  print(summary(object$Dispersion_field_normalized))
}
