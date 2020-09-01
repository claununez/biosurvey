context("Grid from polygons")

test_that("SPDF grid was produced successfully", {
  gfp <- grid_from_region(mx, 200)
  colnam <- colnames(gfp@data)
  gfp_crs <- gfp@proj4string@projargs

  testthat::expect_s4_class(gfp, class = "SpatialPolygonsDataFrame")
  testthat::expect_length(gfp, 99)
  testthat::expect_equal(colnam, c("ID", "Longitude", "Latitude"))
  testthat::expect_equal(gfp_crs, "+proj=longlat +datum=WGS84 +no_defs")
})


test_that("Errors and messages gfr", {
  testthat::expect_error(grid_from_region(mx))
  testthat::expect_error(grid_from_region(cell_size = 100))
  testthat::expect_error(grid_from_region(mx, 5000))
  testthat::expect_error(grid_from_region(mx, c(1000, 1000, 1000)))
  testthat::expect_message(grid_from_region(mx, 1000, complete_cover = FALSE))
})
