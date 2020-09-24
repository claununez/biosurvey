context("Preparation of master matrix")

test_that("Correct creation of master_matrix", {
  variables <- raster::stack(system.file("extdata/variables.tif",
                                         package = "biosurvey"))
  names(variables) <- c("Mean_temperature", "Max_temperature",
                        "Min_temperature", "Annual_precipitation",
                        "Prec_wettest_month", "Prec_driest_month")
  m_matrix <- prepare_master_matrix(region = mx, variables = variables)
  m_matrix1 <- prepare_master_matrix(region = mx, variables = variables,
                                     do_pca = T, scale = T,
                                     variables_in_matrix = c("Mean_temperature",
                                                             "Annual_precipitation"))

  cnam <- names(m_matrix)
  inams <- colnames(m_matrix$data_matrix)
  inams1 <- colnames(m_matrix1$data_matrix)

  anames <- c("data_matrix", "preselected_sites", "region", "mask",
              "raster_base", "PCA_results" )
  enams <- c("Longitude", "Latitude", "Mean_temperature", "Max_temperature",
             "Min_temperature", "Annual_precipitation", "Prec_wettest_month",
             "Prec_driest_month")
  enams1 <- c("Longitude", "Latitude", "Mean_temperature",
              "Annual_precipitation", "PC1", "PC2")


  testthat::expect_s3_class(m_matrix, "master_matrix")
  testthat::expect_null(m_matrix$mask)
  testthat::expect_null(m_matrix$PCA_results)
  testthat::expect_s3_class(m_matrix1$PCA_results, "prcomp")
  testthat::expect_s4_class(m_matrix$region, "SpatialPolygonsDataFrame")
  testthat::expect_s4_class(m_matrix$raster_base, "SpatialPolygonsDataFrame")
  testthat::expect_length(m_matrix, 6)
  testthat::expect_equal(cnam, anames)
  testthat::expect_equal(inams, enams)
  testthat::expect_equal(inams1, enams1)
})


test_that("Errors and messages master_matrix", {
  variables <- raster::stack(system.file("extdata/variables.tif",
                                         package = "biosurvey"))
  testthat::expect_message(prepare_master_matrix(region = mx,
                                                 variables = variables))
  testthat::expect_error(prepare_master_matrix(region = mx,
                                               variables = variables[[1]]))
  testthat::expect_error(prepare_master_matrix(region = mx))
  testthat::expect_error(prepare_master_matrix(variables = variables))
  testthat::expect_error(prepare_master_matrix(region = mx, variables = mx))
  testthat::expect_error(prepare_master_matrix(region = variables,
                                               variables = variables))
})
#----
