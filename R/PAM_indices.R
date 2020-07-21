#' Biodiversity indices derived from PAM
#'
#' @description Calculates a set of biodiversity indices using values contained
#' in a presence-absence matrix.
#'
#' @param PAM matrix, data.frame, or base_PAM object containing information on
#' species presence and absence for a set of sites. Sites are organized in the
#' rows and species in the columns. See details.
#' @param indices (character) code for indices to be calculated. Basic indices
#' are calculated all the time, other indices need to be specified. Options are:
#' "all", "AB", "BW", "BL", "SCSC", "SCSR", "DF", "CC", "WRN", "SRC", "CMSC",
#' and "CMSR". See details. Default = "all"
#' @param exclude_column (optional) name or numeric index of columns to be
#' excluded.Default = NULL.
#'
#' @return
#' If \code{PAM} is a matrix or data.frame, the result is a list with the results
#' described below (depending on \code{indices}). If \code{PAM} is a base_PAM
#' object, a base_PAM object will be returned and the list described above will
#' be appended to the element PAM_indices in such an element.
#'
#' @details
#' Description of the codes of all indices to be calculated is presented in the
#' table below. Basic indices are calculated in all cases not matter the codes
#' defined in \code{indices}. For further details on the way calculations are
#' performed and the meaning of the indices see Soberon and Cavner (2015)
#' \doi{https://doi.org/10.17161/bi.v10i0.4801}.
#'
#' |Code  |Index                                    |Calculation            |
#' |:-----|----------------------------------------:|----------------------:|
#' |      |Richness                                 |Basic                  |
#' |      |Range                                    |Basic                  |
#' |      |Richness standarized                     |Basic                  |
#' |      |Range standarized                        |Basic                  |
#' |AB    |Additive Beta                            |Needs to be defined    |
#' |BW    |Beta Whittaker                           |Needs to be defined    |
#' |BL    |Beta Legendre                            |Needs to be defined    |
#' |SCSC  |Schluter covariance sites-composition    |Needs to be defined    |
#' |SCSR  |Schluter covariance species-ranges       |Needs to be defined    |
#' |DF    |Dispersal field                          |Needs to be defined    |
#' |SCC   |Shared community composition             |Needs to be defined    |
#' |WRN   |Wright-Reeves nestedness                 |Needs to be defined    |
#' |SRC   |Stone-Roberts Cscore                     |Needs to be defined    |
#' |CMSC  |Covariance matrix sites-composition      |Needs to be defined    |
#' |CMSR  |Covariance matrix species-ranges         |Needs to be defined    |
#'
#' @usage
#' PAM_indices(PAM, indices = "all", exclude_column = NULL)
#'
#' @export
#'
#' @examples
#' # Data
#' data("sp_data", package = "biosurvey")
#'
#' # PAM
#' pam <- PAM_from_table(data = sp_data, ID_column = "ID", species_column = "Species")
#'
#' pam_ind <- PAM_indices(pam, exclude_column = 1)
#' pam_ind[1:3]

PAM_indices <- function(PAM, indices = "all", exclude_column = NULL) {
  # Initial test
  if (missing(PAM)) {
    stop("Argument 'PAM' must be defined.")
  }

  all_in <- c("all", "AB", "BW", "BL", "SCSC", "SCSR", "DF", "SCC", "WRN", "SRC",
              "CMSC", "CMSR")
  if (any(!indices %in% all_in)) {
    stop("One or more elements defined in 'indices' is not valid, check function's help.")
  }

  # More tests and preparing data
  cpam <- class(PAM)[1]
  if (!cpam %in% c("base_PAM", "matrix", "data.frame")) {
    stop("Argument 'PAM' must be of class 'base_PAM' or 'matrix'.")
  } else {
    if (cpam == "base_PAM") {
      bpam <- PAM
      PAM <- as.matrix(bpam$PAM@data[, -(1:3)])
      rownames(PAM) <- bpam$PAM@data[, "ID"]
      exclude_column <- NULL
    }
  }

  if (!is.null(exclude_column)) {
    ## test class of exclude column
    if (class(!exclude_column)[1] %in% c("numeric", "character")) {
      stop("Argument 'exclude_column' must be of class 'numeric' or 'character'.")
    }

    ## process matrix to exclude
    if (is.numeric(exclude_column)) {
      PAM <- PAM[, -exclude_column]
    } else {
      PAM <- PAM[, !colnames(PAM) %in% exclude_column]
    }
  }

  if (cpam == "data.frame") {
    PAM <- as.matrix(PAM)
  }

  # PAM transpose
  tm1 <- t(PAM)

  # PAM properties
  S <- ncol(PAM)
  N <- nrow(PAM)

  # matrices A and Omega
  A <- PAM %*% tm1
  O <- tm1 %*% PAM

  # richness and ranges
  ## richness (spp per cell) and ranges (ncells per sp)
  rich <- diag(A)
  rang <- diag(O)

  ## richness and ranges adjusted to S and N
  richS <- rich / S
  rangN <- rang / N

  ## traces of richness and ranges
  trA <- sum(rich)
  trO <- sum(rang)

  # dispersal field and shared community composition
  ## values of dispersal field
  if (any(indices %in% c("all", "DF"))) {
    d_field <- c(PAM %*% rang)
    names(d_field) <- rownames(PAM)
    ## average
    av_dfield <- mean(d_field)
  } else {
    d_field <- av_dfield <- NA
  }

  ## values shared community composition
  if (any(indices %in% c("all", "SCC"))) {
    sc_comp <- c(tm1 %*% rich)
    names(sc_comp) <- colnames(PAM)
    ## average
    av_sccomp <- mean(sc_comp)
  } else {
    sc_comp <- av_sccomp <- NA
  }

  # other indices
  ## Whittaker's multiplicative beta
  if (any(indices %in% c("all", "BW"))) {
    BW <- (S * N) / trO
  } else {
    BW <- NA
  }

  ## Lande's additive beta
  if (any(indices %in% c("all", "AB"))) {
    BA <- S * (1 - (trO / (S * N)))
  } else {
  }

  ## Legendre's beta
  if (any(indices %in% c("all", "BL"))) {
    BL <- trO - sum(d_field)
  } else {
    BL <- NA
  }

  ## Matrix of covariance of composition of sites
  if (any(indices %in% c("all", "CMSC"))) {
    CS_cov <- (A / S) - (richS %*% t(richS))
    ## Mean
    Ccov_mean <- (d_field / (N * S)) - (BW^-1 * richS)
  } else {
    CS_cov <- Ccov_mean <- NA
  }

  ## Matrix of covariance of ranges of species
  if (any(indices %in% c("all", "CMSR"))) {
    RS_cov <- (O / S) - (rangN %*% t(rangN))
    ## Mean
    Rcov_mean <- (sc_comp / (N * S)) - (BW^-1 * rangN)
  } else {
    RS_cov <- Rcov_mean <- NA
  }

  ## Schluter sites-composition covariance
  if (any(indices %in% c("all", "SCSC"))) {
    VCS_cov <- sum(CS_cov) / sum(diag(CS_cov))
  } else {
    VCS_cov <- NA
  }

  ## Schluter species-ranges covariance
  if (any(indices %in% c("all", "SCSR"))) {
    VRS_cov <- sum(RS_cov) / sum(diag(RS_cov))
  } else {
    VRS_cov <- NA
  }

  ## Wright & Reeves' nestedness
  if (any(indices %in% c("all", "WRN"))) {
    Nc <- (sum(d_field) - ((N * S) / BW)) / 2
  } else {
    Nc <- NA
  }

  ## Stone & Roberts Cscore (Hadamard product (x) = element wise multiplication)
  if (any(indices %in% c("all", "SRC"))) {
    Cs <- (sum(O * O) - (N * av_dfield)) / 2
  } else {
    Cs <- NA
  }

  # returning results
  tab_in <- data.frame(Value = c(N, S, av_dfield, av_sccomp, BA, BW, BL, VCS_cov,
                                 VRS_cov, Nc, Cs),
                       row.names = c("Sites_Cells", "Species", "Av_dispersal_field",
                                     "Av_shared_community_composition",
                                     "Additive_Beta", "Beta_Whittaker",
                                     "Beta_Legendre", "Schluter_cov_sites_composition",
                                     "Schluter_cov_species_ranges",
                                     "Wright_Reeves_nestedness",
                                     "Stone_Roberts_Cscore"))

  # if base_PAM
  if (cpam == "base_PAM") {
    bpam$PAM_indices <- list(One_value_indices = tab_in, Richness = rich,
                              Range = rang, Richness_standardized = richS,
                              Range_standardized = rangN, Dispersal_field = d_field,
                              Shared_community_composition = sc_comp,
                              Cov_mat_sites_composition = CS_cov,
                              Cov_mat_species_ranges = RS_cov)
    return(bpam)
  } else {
    # if matrix
    return(list(One_value_indices = tab_in, Richness = rich, Range = rang,
                Richness_standardized = richS, Range_standardized = rangN,
                Dispersal_field = d_field, Shared_comm_composition = sc_comp,
                Cov_mat_sites_composition = CS_cov, Cov_mat_species_ranges = RS_cov))
  }
}
