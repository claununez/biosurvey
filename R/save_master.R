#' Functions to save and read biosurvey objects
#'
#' @param master_matrix object derived from function prepare_master_matrix.
#' @param file_name path (with object name) for saving or reading the master matrix.
#' @param verbose whether or not to print messages about the process.
#' Default = TRUE
#'
#' @return if verbose = TRUE, returns a message indicating the path where the object was saved.
#' @rdname save_master
#' @export
#'
save_master_matrix <- function(master_matrix, file_name, verbose = TRUE) {
  if(class(master_matrix) != "master_matrix") {
    stop("Argument 'master_matrix' must be of class 'master_matrix' or 'master_selection'") }

  #Wrap spatial objects
  if(!is.null(master_matrix$region)) {
  master_matrix$region <- terra::wrap(master_matrix$region) }
  if(!is.null(master_matrix$mask)) {
  master_matrix$mask <- terra::wrap(master_matrix$mask) }
  if(!is.null(master_matrix$raster_base)) {
  master_matrix$raster_base <- terra::wrap(master_matrix$raster_base) }

  saveRDS(master_matrix, paste0(file_name, ".rds"))
  if(verbose) {
  message("Master_matrix saved in ", file_name, ".rds") }

    } #End of function


#' @return a master matrix object with unwrapped spatial objects.
#' @rdname save_master
#' @examples
#' mm <- read_master_matrix(file_name = "data/m_matrix.rds")
read_master_matrix <- function(file_name) {
  master_matrix <- readRDS(file_name)
  #Unwrap spatial objects
  if(!is.null(master_matrix$region)) {
    master_matrix$region <- terra::unwrap(master_matrix$region) }
  if(!is.null(master_matrix$mask)) {
    master_matrix$mask <- terra::unwrap(master_matrix$mask) }
  if(!is.null(master_matrix$raster_base)) {
    master_matrix$raster_base <- terra::unwrap(master_matrix$raster_base) }
  return(master_matrix)
}
