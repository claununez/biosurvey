context("Preparation of mask based on preselected")


test_that("Correct creation of mask", {
  data("m_matrix_pre", package = "biosurvey")
  pdm <- preselected_dist_mask(master = m_matrix_pre, expected_points = 20,
                               space = "G")

  cnam <- names(pdm)
  anames <- c("distance", "mask")
  clth <- class(pdm)

  testthat::expect_equal(clth, "list")
  testthat::expect_s4_class(pdm$mask, "SpatialPolygonsDataFrame")
  testthat::expect_length(pdm$distance, 1)
  testthat::expect_equal(cnam, anames)
})

test_that("Errors and messages preselected_dist_mask", {
  data("m_matrix_pre", package = "biosurvey")

  testthat::expect_error(preselected_dist_mask())
  testthat::expect_error(preselected_dist_mask(master = m_matrix_pre))
  testthat::expect_error(preselected_dist_mask(master = m_matrix_pre,
                                               space = "G"))
  testthat::expect_error(preselected_dist_mask(master = m_matrix_pre,
                                               expected_points = 20))
  testthat::expect_error(preselected_dist_mask(master = m_matrix_pre,
                                               expected_points = 20,
                                               space = "E"))
})
