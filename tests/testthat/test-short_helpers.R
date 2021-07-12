#----
context("Find raster extension")

test_that("Errors match_rformat", {
  match_rformat("GTiff")

  testthat::expect_error(match_rformat())
})
#----


#----
context("Project spatial points")

test_that("Errors wgs84_2aed_laea", {
  data("sp_occurrences", package = "biosurvey")
  sp_occ <- wgs84_2aed_laea(sp_occurrences, longitude = "longitude",
                            latitude = "latitude", which = "EA")

  testthat::expect_error(wgs84_2aed_laea())
  testthat::expect_error(wgs84_2aed_laea(longitude = "longitude",
                                         latitude = "latitude", which = "EA"))
  testthat::expect_error(wgs84_2aed_laea(sp_occurrences, latitude = "latitude",
                                         which = "EA"))
  testthat::expect_error(wgs84_2aed_laea(sp_occurrences, longitude = "longitude",
                                         which = "EA"))
})
#----
