


PAM_indices <- function(PAM, indices = "all", exclude_column = NULL) {
  # Initial test
  if (missing(PAM)) {
    stop("Argument 'PAM' must be defined.")
  }

  all_in <- c("all", "AB", "BW", "BL", "SCSC", "SCSR", "DF", "CC", "WRN", "SRC",
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
  if (any(indices %in% c("all", "SC"))) {
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
                              Range = rang, Richness_standard = richS,
                              Range_standard = rangN, Dispersal_field = d_field,
                              Shared_community_composition = sc_comp,
                              Cov_mat_sites_composition = CS_cov,
                              Cov_mat_species_ranges = RS_cov)
    return(bpam)
  } else {
    # if matrix
    return(list(One_value_indices = tab_in, Richness = rich, Range = rang,
                Richness_standard = richS, Range_standard = rangN,
                Dispersal_field = d_field, Shared_comm_composition = sc_comp,
                Cov_mat_sites_composition = CS_cov, Cov_mat_species_ranges = RS_cov))
  }
}
