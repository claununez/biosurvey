
new_master_matrix <- function(data_matrix, preselected_sites = NULL, region,
                              mask = NULL, raster_base, PCA_results = NULL) {
  stopifnot(is.data.frame(data_matrix))
  preclass <- class(preselected_sites)[1]
  reclass <- class(region)[1]
  mclass <- class(mask)[1]
  raclass <- class(raster_base)[1]
  pcaclass <- class(PCA_results)[1]

  if (!preclass %in% c("data.frame", "NULL")){
    stop("'preselected_sites' must be of class 'data.frame', or NULL.")
  }
  if (!reclass %in% c("SpatialPolygons", "SpatialPolygonsDataFrame")) {
    stop("'region' must be of class 'SpatialPolygons' or 'SpatialPolygonsDataFrame'.")
  }
  if (!mclass %in% c("SpatialPolygons", "SpatialPolygonsDataFrame", "NULL")) {
    stop("'mask' must be of class 'SpatialPolygons' or 'SpatialPolygonsDataFrame', or NULL.")
  }
  if (!raclass %in% c("RasterLayer")) {
    stop("'raster_base' must be of class 'RasterLayer'.")
  }
  if (!pcaclass %in% c("prcomp", "NULL")) {
    stop("'PCA_results' must be of class 'prcomp', or NULL.")
  }

  val <- list(data_matrix = data_matrix, preselected_sites = preselected_sites,
              region = region, mask = mask, raster_base = raster_base,
              PCA_results = PCA_results)
  class(val) <- "master_matrix"
  return(val)
}


new_master_selection <- function(data_matrix, preselected_sites = NULL, region,
                                 mask = NULL, raster_base, PCA_results = NULL,
                                 selected_sites_random = NULL, selected_sites_G = NULL,
                                 selected_sites_E = NULL, selected_sites_EG = NULL) {

  randclass <- class(selected_sites_random)[1]
  gclass <- class(selected_sites_G)[1]
  eclass <- class(selected_sites_E)[1]
  egclass <- class(selected_sites_EG)[1]

  if (!randclass %in% c("list", "NULL")) {
    stop("'selected_sites_random' must be of class 'list', or NULL.")
  }
  if (!gclass %in% c("list", "NULL")) {
    stop("'selected_sites_G' must be of class 'list', or NULL.")
  }
  if (!eclass %in% c("list", "NULL")) {
    stop("'selected_sites_E' must be of class 'list', or NULL.")
  }
  if (!egclass %in% c("list", "NULL")) {
    stop("'selected_sites_EG' must be of class 'list', or NULL.")
  }

  val <- new_master_matrix(data_matrix = data_matrix,
                           preselected_sites = preselected_sites, region = region,
                           mask = mask, raster_base = raster_base,
                           PCA_results = PCA_results)
  val <- c(val, list(selected_sites_random = selected_sites_random,
                     selected_sites_G = selected_sites_G,
                     selected_sites_E = selected_sites_E,
                     selected_sites_EG = selected_sites_EG))
  class(val) <- c("master_selection", "master_matrix")
  return(val)
}
