context("Preparation of base PAM")

test_that("Correct creation of base_PAM", {
  b_pam <- prepare_base_PAM(data = species_data, master_matrix = m_matrix,
                            cell_size = 200)

  cnam <- names(b_pam)
  inams <- names(b_pam$PAM_indices)

  enams <- c("One_value_indices", "Richness", "Range", "Richness_normalized",
             "Range_normalized", "Dispersion_field",
             "Shared_community_composition", "Mean_composition_covariance",
             "Mean_range_covariance", "Cov_mat_sites_composition",
             "Cov_mat_species_ranges")

  testthat::expect_s3_class(b_pam, "base_PAM")
  testthat::expect_s4_class(b_pam$PAM, "SpatialPolygonsDataFrame")
  testthat::expect_length(b_pam, 2)
  testthat::expect_length(b_pam$PAM_indices, 11)
  testthat::expect_equal(cnam, c("PAM", "PAM_indices"))
  testthat::expect_equal(inams, enams)
})


test_that("Errors and messages base_PAM", {
  testthat::expect_message(prepare_base_PAM(data = species_data,
                                            master_matrix = m_matrix,
                                            cell_size = 500))
  testthat::expect_error(prepare_base_PAM(data = species_data))
  testthat::expect_error(prepare_base_PAM(master_matrix = m_matrix))
  testthat::expect_error(prepare_base_PAM(cell_size = 1000))
})
#----
