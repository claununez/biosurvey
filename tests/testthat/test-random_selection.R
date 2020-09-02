context("Random selection of sites")

test_that("Correct random master_selection", {
  selection <- random_selection(m_matrix_pre, n_sites = 10, n_samplings = 5)
  cnam <- names(selection)
  nsel <- nrow(selection$selected_sites_random[[1]])
  cls <- class(selection$selected_sites_random)[1]

  anames <- c("data_matrix", "preselected_sites", "region", "mask",
              "raster_base", "PCA_results", "selected_sites_random",
              "selected_sites_G", "selected_sites_E", "selected_sites_EG" )

  testthat::expect_s3_class(selection, "master_selection")
  testthat::expect_null(selection$selected_sites_G)
  testthat::expect_null(selection$selected_sites_E)
  testthat::expect_null(selection$selected_sites_EG)
  testthat::expect_s3_class(selection$selected_sites_random[[1]], "data.frame")
  testthat::expect_length(selection, 10)
  testthat::expect_length(selection$selected_sites_random, 5)
  testthat::expect_equal(cls, "list")
  testthat::expect_equal(cnam, anames)
  testthat::expect_equal(nsel, 10)
})


test_that("Errors and messages random selection", {
  testthat::expect_message(random_selection(m_matrix_pre, n_sites = 10))
  testthat::expect_error(random_selection(1:100, n_sites = 10))
  testthat::expect_error(random_selection(m_matrix))
  testthat::expect_error(random_selection())
  testthat::expect_error(random_selection(n_sites = 10))
})
#----
