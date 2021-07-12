context("Preparation of diversity-range diagram-data")

test_that("Correct creation of PAM_CS", {
  data("b_pam", package = "biosurvey")
  pcs <- prepare_PAM_CS(PAM = b_pam)


  cnam <- names(pcs$PAM_indices$CS_diagram)
  dm <- dim(pcs$PAM_indices$CS_diagram$Randomized_DF)
  clth <- class(pcs$PAM_indices$CS_diagram$Theoretical_boundaries)
  clrd <- class(pcs$PAM_indices$CS_diagram$Randomized_DF)[1]

  anames <- c("Species", "Sites_cells", "Beta_W", "Spearman_cor",
              "Theoretical_boundaries", "Richness_normalized",
              "Dispersion_field_normalized", "S_significance_id",
              "Randomized_DF")

  testthat::expect_s3_class(pcs, "base_PAM")
  testthat::expect_s3_class(pcs$PAM_indices$CS_diagram, "PAM_CS")
  testthat::expect_equal(clrd, "matrix")
  testthat::expect_equal(clth, "list")
  testthat::expect_equal(dm, c(1, 1))
  testthat::expect_length(pcs$PAM_indices$CS_diagram, 9)
  testthat::expect_length(pcs$PAM_indices$CS_diagram$Theoretical_boundaries, 2)
  testthat::expect_equal(cnam, anames)
})

test_that("Errors and messages PAM_CS", {
  data("b_pam", package = "biosurvey")

  testthat::expect_error(prepare_PAM_CS())
  testthat::expect_error(prepare_PAM_CS(PAM = b_pam$PAM_indices))
})
