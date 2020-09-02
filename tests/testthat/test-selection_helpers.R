#----
context("Point thinning")

test_that("Correct point_thinning", {
  data1 <- m_matrix$data_matrix
  data1 <- data1[sample(nrow(data1), 500), ]
  thin <- point_thinning(data1, x_column = "Longitude", y_column = "Latitude",
                         thinning_distance = 200, space = "G",
                         max_n_samples = 1, replicates = 1, set_seed = 1)


  dtda <- dim(thin[[1]])

  testthat::expect_true(dtda[1] < 500)
  testthat::expect_equal(dtda[2], ncol(data1))
  testthat::expect_length(thin, 1)
})


test_that("Errors point_thinning", {
  testthat::expect_error(point_thinning())
  testthat::expect_error(point_thinning(data1))
  testthat::expect_error(point_thinning(data1, x_column = "Longitude"))
  testthat::expect_error(point_thinning(data1, y_column = "Latitude"))
})
#----

#----
context("Obtaining closest to centroid points")

test_that("Correct closest_to_centroid", {
  data1 <- m_matrix$data_matrix
  data1 <- data1[sample(nrow(data1), 1000), ]
  centroid <- closest_to_centroid(data1, x_column = "Longitude",
                                  y_column = "Latitude", space = "G",
                                  n = 1)
  centroid1 <- closest_to_centroid(data1, x_column = "Longitude",
                                   y_column = "Latitude", space = "G",
                                   n = 4)
  dtda <- dim(centroid)
  dtda1 <- dim(centroid1)
  cnam <- colnames(centroid)

  testthat::expect_equal(dtda, c(1, 11))
  testthat::expect_equal(dtda1, c(4, 11))
  testthat::expect_equal(cnam, c(colnames(data1), "id_column"))
})


test_that("Errors closest_to_centroid", {
  testthat::expect_error(closest_to_centroid())
  testthat::expect_error(closest_to_centroid(data1))
  testthat::expect_error(closest_to_centroid(data1, x_column = "Longitude"))
  testthat::expect_error(closest_to_centroid(data1, y_column = "Latitude"))
})
#----


#----
context("Filter sites by distance")

test_that("Correct distance_filter", {
  slist <- m_selection$selected_sites_random
  max_sites <- distance_filter(slist, median_distance_filter = "max")

  cls <- class(max_sites)
  dtda <- dim(max_sites$selection_1)

  testthat::expect_equal(cls, "list")
  testthat::expect_length(max_sites, 1)
  testthat::expect_equal(dtda, c(20, 11))
})


test_that("Errors distance_filter", {
  testthat::expect_error(distance_filter())
  testthat::expect_error(distance_filter(median_distance_filter = "min"))
})
#----
