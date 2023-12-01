#' Functions to save and read master objects
#'
#' @param PAM object of class base_PAM or PAM_subset derived from the functions
#' \code{\link{prepare_base_PAM}} or \code{\link{subset_PAM}}.
#' @param file_name (character) name for the file to save or read the PAM
#' object (includes extension ".rds").
#' @param verbose whether or not to print messages about the process.
#' Default = TRUE.
#'
#' @return
#' If \code{verbose} = TRUE, a message indicating the path where the object
#' was saved.
#'
#' @rdname save_PAM
#' @export

save_PAM <- function(PAM, file_name, verbose = TRUE) {
  if (missing(PAM)) {
    stop("Argument 'PAM' must be defined.")
  }
  if (missing(file_name)) {
    stop("Argument 'file_name' must be defined.")
  }

  if (!class(PAM)[1] %in% c("base_PAM", "PAM_subset")) {
    stop("Argument 'PAM' must be of class 'base_PAM' or 'PAM_subset'.")
  }

  #Wrap spatial objects
  PAM$PAM <- terra::wrap(PAM$PAM)

  saveRDS(PAM, file = file_name)

  if (verbose) {
    message("PAM saved in ", file_name)
  }
}


#' @return an object of class base_PAM with its spatial object unwrapped.
#' @rdname save_PAM
#' @export
#' @examples
#' b_pam <- read_PAM(file_name = "data/b_pam.rds")

read_PAM <- function(file_name) {
  if (missing(file_name)) {
    stop("Argument 'file_name' must be defined.")
  }
  PAM <- readRDS(file_name)

  # Unwrap spatial objects
  PAM$PAM <- terra::unwrap(PAM$PAM)

  return(PAM)
}
