context("PAM biodiversity indices")

test_that("Correct indices list", {
  pam_ind <- PAM_indices(b_pam)

  cnam <- names(pam_ind)
  inams <- names(pam_ind$PAM_indices)

  enams <- c("One_value_indices", "Richness", "Range", "Richness_normalized",
             "Range_normalized", "Dispersion_field", "Shared_community_composition",
             "Mean_composition_covariance", "Mean_range_covariance",
             "Cov_mat_sites_composition", "Cov_mat_species_ranges")

  testthat::expect_s3_class(pam_ind, "base_PAM")
  testthat::expect_length(pam_ind, 2)
  testthat::expect_length(pam_ind$PAM_indices, 11)
  testthat::expect_equal(cnam, c("PAM", "PAM_indices"))
  testthat::expect_equal(inams, enams)
})


test_that("Errors PAM_indices", {
  testthat::expect_error(PAM_indices())
  testthat::expect_error(PAM_indices(1:10))
  testthat::expect_error(PAM_indices(b_pam, indices = "XRF"))
})
#----
