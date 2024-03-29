#' Example of a master_matrix object with no preselected sites
#'
#' An S3 object of class master_matrix.
#' See function \code{\link{prepare_master_matrix}}.
#'
#' @format A list of 6 elements:
#' \describe{
#'   \item{data_matrix}{data.frame with 6276 rows and 10 columns}
#'   \item{preselected_sites}{NULL}
#'   \item{region}{object of class SpatVector}
#'   \item{mask}{NULL}
#'   \item{raster_base}{object of class SpatRaster}
#'   \item{PCA_results}{list of length 4}
#' }
#'
#' @name m_matrix
#'
#' @examples
#' m_matrix <- read_master(system.file("extdata/m_matrix.rds",
#'                                     package = "biosurvey"))
#'
#' print(m_matrix)
NULL


#' Example of a master_matrix object containing preselected sites
#'
#' A S3 object of class master_matrix.
#' See function \code{\link{prepare_master_matrix}}.
#'
#' @format A list of 6 elements:
#' \describe{
#'   \item{data_matrix}{data.frame with 6276 rows and 10 columns}
#'   \item{preselected_sites}{data.frame with 5 rows and 11 columns}
#'   \item{region}{object of class SpatVector}
#'   \item{mask}{NULL}
#'   \item{raster_base}{object of class SpatRaster}
#'   \item{PCA_results}{list of length 4}
#' }
#'
#' @name m_matrix_pre
#'
#' @examples
#' m_matrix_pre <- read_master(system.file("extdata/m_matrix_pre.rds",
#'                                         package = "biosurvey"))
#'
#' print(m_matrix_pre)
NULL


#' Example of a data.frame of preselected sites
#'
#' @description A data.frame with 5 rows and three columns: "Site", "Longitude",
#' and "Latitude".
#'
#' @format data.frame:
#' \describe{
#'   \item{Site}{name of preselected sites}
#'   \item{Longitude}{x coordinates}
#'   \item{Latitude}{y coordinates}
#' }
#'
#' @examples
#' data("preselected", package = "biosurvey")
#'
#' print(preselected)
"preselected"


#' Example of spatial polygon for a region of interest
#'
#' An object of class SpatVector.
#'
#' @format SpatVector
#'
#' @name mx
#'
#' @examples
#' mx <- terra::vect(system.file("extdata/mx.gpkg", package = "biosurvey"))
#'
#' mx
NULL


#' Example of a data.frame of species' found in distinct positions
#'
#' @description A data.frame with 590 rows and two columns: "ID" and "Species".
#' It is the output of the function \code{\link{spdf_2data}}.
#'
#' @format data.frame:
#' \describe{
#'   \item{ID}{identifier of position}
#'   \item{Species}{different species in the table}
#' }
#'
#' @examples
#' data("sp_data", package = "biosurvey")
#'
#' summary(sp_data)
"sp_data"


#' Example of species ranges as a SpatVector
#'
#' An object of class SpatVector.
#'
#' @format SpatVector
#'
#' @name species_data
#'
#' @examples
#' species_data <- terra::vect(system.file("extdata/species_data.gpkg",
#'                                         package = "biosurvey"))
#'
#' species_data
NULL


#' Example of raster layers of suitable and unsuitable conditions for species
#'
#' An object of class SpatRaster containing information about suitable and
#' unsuitable conditions for five species.
#'
#' @format SpatRaster with 109 rows, 182 columns, 19838 cells, and 5 layers:
#' \describe{
#'   \item{Layer}{suitable (1) and unsuitable (0) conditions}
#' }
#'
#' @examples
#' sp_layers <- terra::rast(system.file("extdata/sp_layers.tif",
#'                                      package = "biosurvey"))
#'
#' sp_layers
#' @name sp_layers
NULL


#' Example of variables to be used for preparing a master matrix
#'
#' A dataset containing raster variables for an area that is relevant for used
#' in examples included in the package \code{\link{biosurvey}}.
#'
#' @format A SpatRaster with 109 rows, 190 columns, 20710 cells, and 6 layers:
#' \describe{
#'   \item{Temperature}{temperature, in degrees Celsius times 10.}
#'   \item{Precipitation}{precipitation, in millimeters.}
#' }
#'
#' @source \url{https://www.worldclim.org}
#'
#' @examples
#' variables <- terra::rast(system.file("extdata/variables.tif",
#'                                      package = "biosurvey"))
#'
#' terra::plot(variables[[1]])
#' @name variables
NULL


#' Occurrence records for the species Parides gundlachianus
#'
#' @description Dataset containing geographic coordinates of a Cuban
#' butterfly.
#'
#' @format A data.frame with 19 rows and 3 columns.
#' \describe{
#'   \item{name}{character, species scientific name.}
#'   \item{longitude}{numeric, longitude values.}
#'   \item{latitude}{numeric, latitude values.}
#' }
#' @source \url{https://www.gbif.org/}
#'
#' @examples
#' data("sp_occurrences", package = "biosurvey")
#'
#' head(sp_occurrences)
"sp_occurrences"


#' A list of vectors of distances
#'
#' @description A list of six vectors of point distances.
#'
#' @format A list.
#' \describe{
#'   \item{vector}{numeric, values of distances (six elements)}
#' }
#'
#' @examples
#' data("dist_list", package = "biosurvey")
"dist_list"


#' Example of a master_selection object from using functions for selecting sites
#'
#' An S3 object of class master_selection. See functions
#' \code{\link{uniformE_selection}}, \code{\link{uniformG_selection}},
#' \code{\link{random_selection}}, or \code{\link{EG_selection}}.
#'
#' @format A list of 10 elements:
#' \describe{
#'   \item{data_matrix}{data.frame with 6276 rows and 10 columns}
#'   \item{preselected_sites}{NULL}
#'   \item{region}{object of class SpatVector}
#'   \item{mask}{NULL}
#'   \item{raster_base}{object of class SpatRaster}
#'   \item{PCA_results}{list of length 4}
#'   \item{selected_sites_random}{list with one data.frame}
#'   \item{selected_sites_G}{list with one data.frame}
#'   \item{selected_sites_E}{list with one data.frame}
#'   \item{selected_sites_EG}{NULL}
#' }
#'
#' @name m_selection
#'
#' @examples
#' m_selection <- read_master(system.file("extdata/m_selection.rds",
#'                                        package = "biosurvey"))
#'
#' print(m_selection)
NULL


#' Example of object obtained from using the function base_PAM
#'
#' An S3 object of class base_PAM. See functions \code{\link{prepare_base_PAM}}.
#'
#' @format A list of 2 elements:
#' \describe{
#'   \item{PAM}{SpatVector with 306 features}
#'   \item{PAM_indices}{a list of 11 elements}
#' }
#'
#' @name b_pam
#'
#' @examples
#' b_pam <- read_PAM(system.file("extdata/b_pam.rds",
#'                               package = "biosurvey"))
#'
#' print(b_pam)
NULL

