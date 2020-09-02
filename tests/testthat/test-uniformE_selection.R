context("Uniform selection of sites in E")

test_that("Correct E master_selection", {
  m_blocks <- make_blocks(m_matrix_pre, variable_1 = "PC1",
                          variable_2 = "PC2", n_cols = 10, n_rows = 10,
                          block_type = "equal_area")
  selection <- uniformE_selection(m_blocks, variable_1 = "PC1", variable_2 = "PC2",
                                  selection_from = "block_centroids",
                                  expected_points = 15, max_n_samplings = 1,
                                  replicates = 1)
  cnam <- names(selection)
  nsel <- nrow(selection$selected_sites_E[[1]])
  cls <- class(selection$selected_sites_E)[1]

  anames <- c("data_matrix", "preselected_sites", "region", "mask",
              "raster_base", "PCA_results", "selected_sites_random",
              "selected_sites_G", "selected_sites_E", "selected_sites_EG" )

  testthat::expect_s3_class(selection, "master_selection")
  testthat::expect_null(selection$selected_sites_random)
  testthat::expect_null(selection$selected_sites_G)
  testthat::expect_null(selection$selected_sites_EG)
  testthat::expect_s3_class(selection$selected_sites_E[[1]], "data.frame")
  testthat::expect_length(selection, 10)
  testthat::expect_length(selection$selected_sites_E, 1)
  testthat::expect_equal(cls, "list")
  testthat::expect_equal(cnam, anames)
  testthat::expect_equal(nsel, 15)
})


test_that("Errors and messages E selection", {
  m_blocks <- make_blocks(m_matrix_pre, variable_1 = "PC1",
                          variable_2 = "PC2", n_cols = 10, n_rows = 10,
                          block_type = "equal_area")
  testthat::expect_message(uniformE_selection(m_blocks, variable_1 = "PC1",
                                              variable_2 = "PC2",
                                              selection_from = "block_centroids",
                                              expected_points = 15,
                                              max_n_samplings = 1,
                                              replicates = 1))
  testthat::expect_error(uniformE_selection(m_matrix_pre, variable_1 = "PC1",
                                              variable_2 = "PC2",
                                              selection_from = "block_centroids",
                                              expected_points = 15,
                                              max_n_samplings = 1,
                                              replicates = 1))
  testthat::expect_error(uniformE_selection(1:100, variable_1 = "PC1",
                                            variable_2 = "PC2",
                                            selection_from = "block_centroids",
                                            expected_points = 15,
                                            max_n_samplings = 1,
                                            replicates = 1))
  testthat::expect_error(uniformE_selection(m_matrix_pre,
                                            selection_from = "block_centroids",
                                            expected_points = 15,
                                            max_n_samplings = 1,
                                            replicates = 1))
  testthat::expect_error(uniformE_selection(m_blocks))
  testthat::expect_error(uniformE_selection())
  testthat::expect_error(uniformE_selection(expected_points = 10))
})
#----
