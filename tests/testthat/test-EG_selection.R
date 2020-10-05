context("Selection of sites considering EG")

test_that("Correct EG master_selection", {
  m_blocks <- make_blocks(m_matrix_pre, variable_1 = "PC1",
                          variable_2 = "PC2", n_cols = 10, n_rows = 10,
                          block_type = "equal_area")
  selection <- EG_selection(master = m_blocks, variable_1 = "PC1",
                            variable_2 = "PC2", n_blocks = 6, replicates = 1,
                            max_n_samplings = 1, select_point = "E_centroid",
                            cluster_method = "hierarchical",
                            sample_for_distance = 100)
  cnam <- names(selection)
  nsel <- nrow(selection$selected_sites_EG[[1]])
  cls <- class(selection$selected_sites_EG)[1]

  anames <- c("data_matrix", "preselected_sites", "region", "mask",
              "raster_base", "PCA_results", "selected_sites_random",
              "selected_sites_G", "selected_sites_E", "selected_sites_EG" )

  testthat::expect_s3_class(selection, "master_selection")
  testthat::expect_null(selection$selected_sites_random)
  testthat::expect_null(selection$selected_sites_G)
  testthat::expect_null(selection$selected_sites_E)
  testthat::expect_s3_class(selection$selected_sites_EG[[1]], "data.frame")
  testthat::expect_length(selection, 10)
  testthat::expect_length(selection$selected_sites_EG, 1)
  testthat::expect_equal(cls, "list")
  testthat::expect_equal(cnam, anames)
  testthat::expect_equal(nsel, 7)
})


test_that("Errors and messages EG selection", {
  m_blocks <- make_blocks(m_matrix_pre, variable_1 = "PC1",
                          variable_2 = "PC2", n_cols = 10, n_rows = 10,
                          block_type = "equal_area")
  selection <- EG_selection(master = m_blocks, variable_1 = "PC1",
                            variable_2 = "PC2", n_blocks = 6, replicates = 1,
                            max_n_samplings = 1, select_point = "E_centroid",
                            cluster_method = "hierarchical",
                            sample_for_distance = 100)
  testthat::expect_message(EG_selection(master = m_blocks, variable_1 = "PC1",
                                        variable_2 = "PC2", n_blocks = 6,
                                        replicates = 1, max_n_samplings = 1,
                                        select_point = "E_centroid",
                                        cluster_method = "hierarchical",
                                        sample_for_distance = 100))
  testthat::expect_error(EG_selection(master = m_blocks, variable_1 = "PC1",
                                      n_blocks = 6,
                                      replicates = 1, max_n_samplings = 1,
                                      select_point = "E_centroid",
                                      cluster_method = "hierarchical",
                                      sample_for_distance = 100))
  testthat::expect_error(EG_selection(master = m_blocks,
                                      variable_2 = "PC2", n_blocks = 6,
                                      replicates = 1, max_n_samplings = 1,
                                      select_point = "E_centroid",
                                      cluster_method = "hierarchical",
                                      sample_for_distance = 100))
  testthat::expect_error(EG_selection(master = m_blocks, variable_1 = "PC1",
                                      variable_2 = "PC2",
                                      replicates = 1, max_n_samplings = 1,
                                      select_point = "E_centroid",
                                      cluster_method = "hierarchical",
                                      sample_for_distance = 100))
  testthat::expect_error(EG_selection(m_blocks))
  testthat::expect_error(EG_selection())
  testthat::expect_error(EG_selection(n_blocks = 10))
})
#----

