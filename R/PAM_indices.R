all

Sites-Cells
Species
Richness
Ranges
Dispersal_field
Community_composition
Additive_Beta
Beta_Whittaker
Beta_Legendre
Schluter_cov_sites-composition
Schluter_cov_species-ranges
Wright-Reeves_nestedness
Stone_Roberts_Cscore

Richness_standard
Range_standard
Dispersal_field
Sha_comm_composition
Cov_mat_composition_sites
Cov_mat_ranges_species

m1 <- matrix(sample(c(0, 1), 300, replace = T), 20, 15)

PAM_indices <- function(PAM, indices = "all") {
  # Initial test
  if (missing(PAM)) {
    stop("Argument 'PAM' must be defined.")
  }
  
  all_in <- c("all", "AB", "BW", "BL", "SCSC", "SCSR", "DF", "CC", "WRN", "SRC", 
              "CMSC", "CMSR")
  if (any(!indices %in% all_in)) {
    stop("One or more elements defined in 'indices' is not valid, check function's help.")
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
    ## average
    av_dfield <- mean(d_field) 
  } else {
    d_field <- av_dfield <- NA
  }
  
  ## values shared community composition
  if (any(indices %in% c("all", "SC"))) {
    sc_comp <- c(tm1 %*% rich)
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
  tab_in <- data.frame(Index_metric = c("Sites_Cells", "Species", "Av_dispersal_field",
                                        "Av_shared_comm_composition", "Additive_Beta", 
                                        "Beta_Whittaker", "Beta_Legendre", 
                                        "Schluter_cov_sites_composition",
                                        "Schluter_cov_species_ranges",
                                        "Wright_Reeves_nestedness", 
                                        "Stone_Roberts_Cscore"), 
                       Value = c(N, S, av_dfield, av_sccomp, BA, BW, BL, VCS_cov,
                                 VRS_cov, Nc, Cs))
  
  return(list(One_value_indices = tab_in, Richness = rich, Range = rang, 
              Richness_standard = richS, Range_standard = rangN, 
              Dispersal_field = d_field, Shared_comm_composition = sc_comp,
              Cov_mat_sites_composition = CS_cov, Cov_mat_species_ranges = RS_cov))
}

test <- PAM_indices(m1)
test
