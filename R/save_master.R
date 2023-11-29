#' Functions to save and read master objects
#'
#' @param master object derived from the functions
#' \code{\link{prepare_master_matrix}}, \code{\link{random_selection}},
#' \code{\link{uniformG_selection}}, \code{\link{uniformE_selection}}, or
#' \code{uniformEG_selection}.
#' @param file_name (character) name for the file to save or read the master
#' object.
#' @param verbose whether or not to print messages about the process.
#' Default = TRUE.
#'
#' @return
#' If \code{verbose} = TRUE, a message indicating the path where the object
#' was saved.
#'
#' @rdname save_master
#' @export

save_master <- function(master, file_name, verbose = TRUE) {
  if (missing(master)) {
    stop("Argument 'master' must be defined.")
  }
  if (missing(file_name)) {
    stop("Argument 'file_name' must be defined.")
  }

  if (class(master)[1] != "master") {
    stop("Argument 'master' must be of class 'master' or 'master_selection'.")
  }

  #Wrap spatial objects
  master$region <- terra::wrap(master$region)
  master$raster_base <- terra::wrap(master$raster_base)

  if (!is.null(master$mask)) {
    master$mask <- terra::wrap(master$mask)
  }

  saveRDS(master, paste0(file_name, ".rds"))

  if (verbose) {
    message("master saved in ", file_name, ".rds")
  }
}


#' @return a master matrix object with unwrapped spatial objects.
#' @rdname save_master
#' @export
#' @examples
#' mm <- read_master(file_name = "data/m_matrix.rds")

read_master <- function(file_name) {
  if (missing(file_name)) {
    stop("Argument 'file_name' must be defined.")
  }
  master <- readRDS(file_name)

  # Unwrap spatial objects
  master$region <- terra::unwrap(master$region)
  master$raster_base <- terra::unwrap(master$raster_base)

  if(!is.null(master$mask)) {
    master$mask <- terra::unwrap(master$mask)
  }

  return(master)
}
