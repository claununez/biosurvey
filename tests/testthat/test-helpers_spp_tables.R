#----
context("Tables from stack")

test_that("Correct stack_2data", {
  rsp <- raster::stack(system.file("extdata/sp_layers.tif",
                                   package = "biosurvey"))
  names(rsp) <- paste0("Species_", 1:5)
  sp_data <- stack_2data(species_layers = rsp)

  dtda <- dim(sp_data)
  cnam <- colnames(sp_data)
  uval <- unique(sp_data[, 3])

  testthat::expect_equal(dtda, c(11191, 3))
  testthat::expect_equal(cnam, c("Longitude", "Latitude", "Species"))
  testthat::expect_equal(uval, paste0("Species_", 1:5))
})


test_that("Errors stack_2data", {
  testthat::expect_error(stack_2data())
  testthat::expect_error(stack_2data(species_layers = rsp[[1]]))
})
#----

#----
context("Tables from spdf")

test_that("Correct spdf_2data", {
  grid_reg <- grid_from_region(region = mx, cell_size = 200)
  sp_data <- spdf_2data(spdf_object = species_data, spdf_grid = grid_reg)

  dtda <- dim(sp_data)
  cnam <- colnames(sp_data)
  uval <- levels(sp_data[, 2])

  testthat::expect_equal(dtda, c(268, 2))
  testthat::expect_equal(cnam, c("ID", "Species"))
  testthat::expect_equal(uval, paste0("Species_", 1:25))
})


test_that("Errors spdf_2data", {
  testthat::expect_error(spdf_2data())
  testthat::expect_error(spdf_2data(spdf_object = species_data))
  testthat::expect_error(spdf_2data(spdf_grid = grid_reg))
  testthat::expect_error(spdf_2data(spdf_object = species_data,
                                    spdf_grid = as(grid_reg,
                                                   "SpatialPolygons")))
})
#----


#----
context("Tables from rlist")

test_that("Correct rlist_2data", {
  rsp <- raster::stack(system.file("extdata/sp_layers.tif",
                                   package = "biosurvey"))
  names(rsp) <- paste0("Species_", 1:5)
  rlist <- lapply(1:5, function(x) {rsp[[x]]})
  sp_data <- rlist_2data(raster_list = rlist)

  dtda <- dim(sp_data)
  cnam <- colnames(sp_data)
  uval <- unique(sp_data[, 3])

  testthat::expect_equal(dtda, c(11191, 3))
  testthat::expect_equal(cnam, c("Longitude", "Latitude", "Species"))
  testthat::expect_equal(uval, paste0("Species_", 1:5))
})


test_that("Errors rlist_2data", {
  testthat::expect_error(rlist_2data())
  testthat::expect_error(rlist_2data(raster_list = rlist[[1]]))
  testthat::expect_error(rlist_2data(raster_list = rsp))
})
#----
