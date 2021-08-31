context("Selecting E blocks")

test_that("Correct selection of blocks", {
  m_blocks <- make_blocks(m_matrix, variable_1 = "PC1",
                          variable_2 = "PC2", n_cols = 10, n_rows = 10,
                          block_type = "equal_area")
  block_sel <- block_sample(m_blocks, expected_blocks = 10,
                            selection_type = "uniform")

  cnam <- names(block_sel)
  inams <- colnames(block_sel$data_matrix)

  anames <- c("data_matrix", "preselected_sites", "region", "mask",
              "raster_base", "PCA_results" )
  enams <- c("Longitude", "Latitude", "Mean_temperature", "Max_temperature",
             "Min_temperature", "Annual_precipitation", "Prec_wettest_month",
             "Prec_driest_month", "PC1", "PC2", "Block", "Selected_blocks")

  testthat::expect_s3_class(block_sel, "master_matrix")
  testthat::expect_length(block_sel, 6)
  testthat::expect_equal(cnam, anames)
  testthat::expect_equal(inams, enams)
})


test_that("Errors block_sample", {
  m_blocks <- make_blocks(m_matrix, variable_1 = "PC1",
                          variable_2 = "PC2", n_cols = 10, n_rows = 10,
                          block_type = "equal_area")
  testthat::expect_error(block_sample())
  testthat::expect_error(block_sample(m_blocks))
})
#----
