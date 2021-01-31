#' Constructor of S3 objects of class master_matrix
#'
#' @name master_matrix
#' @aliases master_matrix new_master_matrix
#'
#' @param data_matrix a date.frame with information about geographic location of
#' raster cells, initial environmental data, and if available, the first two
#' principal components derived from an analysis done with environmental values.
#' @param preselected_sites data.frame containing sites that must be included
#' in posterior selections of sites for the survey system. Columns must be:
#' "Sites", "Longitude", "Latitude", in that order. Default = NULL.
#' @param region SpatialPolygons* object representing the region of interest.
#' @param mask SpatialPolygons* object used. Default = NULL.
#' @param raster_base a SpatialPolygonsDataFrame representing the grid of the
#' raster layers used, which will be used for plotting purposes.
#' @param PCA_results results of principal component analysis performed with
#' values from raster layers used. Default = NULL.
#'
#' @return
#' An S3 object of class \code{master_matrix}.
#'
#' @export
#'
#' @usage
#' new_master_matrix(data_matrix, preselected_sites = NULL, region,
#'                   mask = NULL, raster_base, PCA_results = NULL)

new_master_matrix <- function(data_matrix, preselected_sites = NULL, region,
                              mask = NULL, raster_base, PCA_results = NULL) {
  if (missing(data_matrix)) {
    stop("Argument 'data_matrix' must be defined")
  }
  if (missing(region)) {
    stop("Argument 'region' must be defined")
  }
  if (missing(raster_base)) {
    stop("Argument 'raster_base' must be defined")
  }
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


#' Constructor of S3 objects of class master_selection
#'
#' @name master_selection
#' @aliases master_selection new_master_selection
#'
#' @param data_matrix a date.frame with information about geographic location of
#' raster cells, initial environmental data, and if available, the first two
#' principal components derived from an analysis done with environmental values.
#' @param preselected_sites data.frame containing sites that must be included
#' in posterior selections of sites for the survey system. Columns must be:
#' "Sites", "Longitude", "Latitude", in that order. Default = NULL.
#' @param region SpatialPolygons* object representing the region of interest.
#' @param mask SpatialPolygons* object used. Default = NULL.
#' @param raster_base a SpatialPolygonsDataFrame representing the grid of the
#' raster layers used, which will be used for plotting purposes.
#' @param PCA_results results of principal component analysis performed with
#' values from raster layers used. Default = NULL.
#' @param selected_sites_random data.frame with the sites selected randomly.
#' Default = NULL.
#' @param selected_sites_G data.frame with the sites selected based on
#' geographic distances. Default = NULL.
#' @param selected_sites_E data.frame with the sites selected based on
#' environmental distances. Default = NULL.
#' @param selected_sites_EG data.frame with the sites selected based on
#' environmental and geographic considerations. Default = NULL.
#'
#' @return
#' An S3 object of class \code{master_selection}.
#'
#' @export
#'
#' @usage
#' new_master_selection(data_matrix, preselected_sites = NULL, region,
#'                      mask = NULL, raster_base, PCA_results = NULL,
#'                      selected_sites_random = NULL, selected_sites_G = NULL,
#'                      selected_sites_E = NULL, selected_sites_EG = NULL)

new_master_selection <- function(data_matrix, preselected_sites = NULL, region,
                                 mask = NULL, raster_base, PCA_results = NULL,
                                 selected_sites_random = NULL,
                                 selected_sites_G = NULL,
                                 selected_sites_E = NULL,
                                 selected_sites_EG = NULL) {

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
                           preselected_sites = preselected_sites,
                           region = region, mask = mask,
                           raster_base = raster_base,
                           PCA_results = PCA_results)
  val <- c(val, list(selected_sites_random = selected_sites_random,
                     selected_sites_G = selected_sites_G,
                     selected_sites_E = selected_sites_E,
                     selected_sites_EG = selected_sites_EG))
  class(val) <- c("master_selection", "master_matrix")
  return(val)
}



#' Constructor of S3 objects of class base_PAM
#'
#' @name base_PAM
#'
#' @param PAM a SpatialPolygonsDataFrame object associated to information about
#' species presence and absence in a goegraphic grid.
#' @param PAM_indices list of indices derived from a PAM. Default = NULL.
#'
#' @export
#'
#' @return
#' An object of class \code{base_PAM}.
#'
#' @usage
#' new_base_PAM(PAM = new("SpatialPolygonsDataFrame"), PAM_indices = NULL)

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


#' Constructor of S3 objects of class PAM_subset
#'
#' @name PAM_subset
#'
#' @param PAM a SpatialPolygonsDataFrame object associated to information about
#' species presence and absence in a goegraphic grid.
#' @param PAM_indices list of indices derived from a PAM. Default = NULL.
#' @param PAM_selected_sites_random subset of \code{PAM} for sites derived from
#' random selection. Default = NULL.
#' @param PAM_selected_sites_G subset of \code{PAM} for sites derived from
#' selection considering geographic distances. Default = NULL.
#' @param PAM_selected_sites_E subset of \code{PAM} for sites derived from
#' selection considering environmental distances. Default = NULL.
#' @param PAM_selected_sites_EG subset of \code{PAM} for sites derived from
#' selection considering environment and geography. Default = NULL.
#'
#' @export
#'
#' @return
#' An object of class \code{PAM_subset}.
#'
#' @usage
#' new_PAM_subset(PAM = new("SpatialPolygonsDataFrame"), PAM_indices = NULL,
#'                PAM_selected_sites_random = NULL, PAM_selected_sites_G = NULL,
#'                PAM_selected_sites_E = NULL, PAM_selected_sites_EG = NULL)

new_PAM_subset <- function(PAM = new("SpatialPolygonsDataFrame"),
                           PAM_indices = NULL, PAM_selected_sites_random = NULL,
                           PAM_selected_sites_G = NULL,
                           PAM_selected_sites_E = NULL,
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



#' Constructor of S3 objects of class PAM_CS
#'
#' @name PAM_CS
#'
#' @param Species (numeric) species name. Default = NA.
#' @param Sites_cells (numeric) number of sites or cells. Default = NA.
#' @param Beta_W (numeric) value of Whittaker's Beta. Default = NA.
#' @param Spearman_cor (numeric) value of Spearman's correlation. Default = NA.
#' @param Theoretical_boundaries list of theoretical boundaries for the values.
#' Default = NA.
#' @param Richness_normalized (numeric) values of normalized richness.
#' Default = NA.
#' @param Dispersion_field_normalized (numeric) values of normalized dispersion
#' field. Default = NA.
#' @param S_significance_id (numeric) values indicating statistical significance
#' of the normalized dispersion field. Default = NA.
#' @param Randomized_DF matrix of values resulted from randomizing matrices.
#' Default = NA.
#'
#' @export
#'
#' @return
#' An object of class \code{PAM_CS}.
#'
#' @usage
#' new_PAM_CS(Species = NA, Sites_cells = NA, Beta_W = NA, Spearman_cor = NA,
#'            Theoretical_boundaries = list(x = NA, y = NA),
#'            Richness_normalized = NA, Dispersion_field_normalized = NA,
#'            S_significance_id = NA, Randomized_DF = matrix())

new_PAM_CS <- function(Species = NA_integer_, Sites_cells = NA_integer_,
                       Beta_W = NA_real_, Spearman_cor = NA_real_,
                       Theoretical_boundaries = list(x = NA_real_, y = NA_real_),
                       Richness_normalized = NA_real_,
                       Dispersion_field_normalized = NA_real_,
                       S_significance_id = NA_integer_,
                       Randomized_DF = matrix()) {

  stopifnot(is.numeric(Species))
  stopifnot(is.numeric(Sites_cells))
  stopifnot(is.numeric(Beta_W))
  stopifnot(is.numeric(Spearman_cor))
  stopifnot(is.list(Theoretical_boundaries))
  stopifnot(is.numeric(Richness_normalized))
  stopifnot(is.numeric(Dispersion_field_normalized))
  stopifnot(is.numeric(S_significance_id))
  stopifnot(is.matrix(Randomized_DF))

  val <- list(Species = Species, Sites_cells = Sites_cells,
              Beta_W = Beta_W, Spearman_cor = Spearman_cor,
              Theoretical_boundaries = Theoretical_boundaries,
              Richness_normalized = Richness_normalized,
              Dispersion_field_normalized = Dispersion_field_normalized,
              S_significance_id = S_significance_id,
              Randomized_DF = Randomized_DF)

  class(val) <- "PAM_CS"
  return(val)
}
