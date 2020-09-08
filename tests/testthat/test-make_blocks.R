context("Creating E blocks")

test_that("Correct creation of blocks", {
  m_blocks <- make_blocks(m_matrix, variable_1 = "PC1",
                          variable_2 = "PC2", n_cols = 10, n_rows = 10,
                          block_type = "equal_area")

  cnam <- names(m_blocks)
  inams <- colnames(m_blocks$data_matrix)

  anames <- c("data_matrix", "preselected_sites", "region", "mask",
              "raster_base", "PCA_results" )
  enams <- c("Longitude", "Latitude", "Mean_temperature", "Max_temperature",
             "Min_temperature", "Annual_precipitation", "Prec_wettest_month",
             "Prec_driest_month", "PC1", "PC2", "Block")

  testthat::expect_s3_class(m_matrix, "master_matrix")
  testthat::expect_length(m_matrix, 6)
  testthat::expect_equal(cnam, anames)
  testthat::expect_equal(inams, enams)
})


test_that("Errors and messages make_blocks", {
  testthat::expect_error(make_blocks())
  testthat::expect_error(make_blocks(master_matrix = m_matrix))
  testthat::expect_error(make_blocks(master_matrix = m_matrix, variable_1 = "PC1"))
  testthat::expect_error(make_blocks(region = variables, variable_1 = "PC1",
                                     variable_2 = "PC2"))
})
#----
