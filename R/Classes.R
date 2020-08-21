
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
  if (!raclass %in% c("SpatialPolygonsDataFrame")) {
    stop("'raster_base' must be of class 'SpatialPolygonsDataFrame'.")
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



new_base_PAM <- function(PAM = new("SpatialPolygonsDataFrame"),
                         PAM_indices = NULL) {
  pclass <- class(PAM)[1]
  piclass <- class(PAM_indices)[1]

  if (!pclass %in% c("SpatialPolygonsDataFrame")) {
    stop("'PAM' must be of class 'list', or NULL.")
  }
  if (!piclass %in% c("list", "NULL")) {
    stop("'PAM_indices' must be of class 'list', or NULL.")
  }

  val <- list(PAM = PAM, PAM_indices = PAM_indices)
  class(val) <- "base_PAM"
  return(val)
}



new_PAM_subset <- function(PAM = new("SpatialPolygonsDataFrame"),
                           PAM_indices = NULL, PAM_selected_sites_random = NULL,
                           PAM_selected_sites_G = NULL, PAM_selected_sites_E = NULL,
                           PAM_selected_sites_EG = NULL) {
  prclass <- class(PAM_selected_sites_random)[1]
  pgclass <- class(PAM_selected_sites_G)[1]
  peclass <- class(PAM_selected_sites_E)[1]
  pegclass <- class(PAM_selected_sites_EG)[1]


  if (!prclass %in% c("list", "NULL")) {
    stop("'PAM_selected_sites_random' must be of class 'list', or NULL.")
  }
  if (!pgclass %in% c("list", "NULL")) {
    stop("'PAM_selected_sites_G' must be of class 'list', or NULL.")
  }
  if (!peclass %in% c("list", "NULL")) {
    stop("'PAM_selected_sites_E' must be of class 'list', or NULL.")
  }
  if (!pegclass %in% c("list", "NULL")) {
    stop("'PAM_selected_sites_EG' must be of class 'list', or NULL.")
  }

  val <- new_base_PAM(PAM, PAM_indices)
  val <- c(val, list(PAM_selected_sites_random = PAM_selected_sites_random,
                     PAM_selected_sites_G = PAM_selected_sites_G,
                     PAM_selected_sites_E = PAM_selected_sites_E,
                     PAM_selected_sites_EG = PAM_selected_sites_EG))
  class(val) <- c("PAM_subset", "base_PAM")
  return(val)
}
