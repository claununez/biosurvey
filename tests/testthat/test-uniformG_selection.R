context("Uniform selection of sites in G")

test_that("Correct G master_selection", {
  selection <- uniformG_selection(m_matrix_pre, expected_points = 20,
                                  max_n_samplings = 1, replicates = 1)
  cnam <- names(selection)
  nsel <- nrow(selection$selected_sites_G[[1]])
  cls <- class(selection$selected_sites_G)[1]

  anames <- c("data_matrix", "preselected_sites", "region", "mask",
              "raster_base", "PCA_results", "selected_sites_random",
              "selected_sites_G", "selected_sites_E", "selected_sites_EG" )

  testthat::expect_s3_class(selection, "master_selection")
  testthat::expect_null(selection$selected_sites_random)
  testthat::expect_null(selection$selected_sites_E)
  testthat::expect_null(selection$selected_sites_EG)
  testthat::expect_s3_class(selection$selected_sites_G[[1]], "data.frame")
  testthat::expect_length(selection, 10)
  testthat::expect_length(selection$selected_sites_G, 1)
  testthat::expect_equal(cls, "list")
  testthat::expect_equal(cnam, anames)
  testthat::expect_equal(nsel, 20)
})


test_that("Errors and messages G selection", {
  testthat::expect_message(uniformG_selection(m_matrix_pre,
                                              expected_points = 20,
                                              max_n_samplings = 1,
                                              replicates = 1))
  testthat::expect_error(uniformG_selection(1:100, expected_points = 20,
                                            max_n_samplings = 1,
                                            replicates = 1))
  testthat::expect_error(uniformG_selection(m_matrix))
  testthat::expect_error(uniformG_selection())
  testthat::expect_error(uniformG_selection(expected_points = 10))
})
#----
