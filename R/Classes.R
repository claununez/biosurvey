new_master_matrix <- function(data_matrix, region, mask, raster_base,
                              preselected_sites, PCA_results) {
  stopifnot(is.data.frame(data_matrix))

  reclass <- class(region)[1]
  mclass <- class(mask)[1]
  raclass <- class(raster_base)[1]
  pcaclass <- class(PCA_results)[1]
  preclass <- class(preselected_sites)[1]

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

  return(structure(list(data_matrix = data_matrix, region = region, mask = mask,
                        raster_base = raster_base,
                        preselected_sites = preselected_sites,
                        PCA_results = PCA_results), class = "master_matrix"))
}
