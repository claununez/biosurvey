
#print.master_matrix <- function(x) {
#  cat("A master_matrix object:")
#
#}


#' Summary of attributes and results
#' @name summary
#' @aliases summary,master_matrix-method
#' @param object object of class master_matrix.
#' @export
#' @return
#' A printed summary.
#' @rdname summary


summary.master_matrix <- function(object) {
  # -----------
  # detecting potential errors
  if (!missing(object)) {
    clo <- class(object)[1]
    if (clo != "master_matrix") {
      stop("Argument 'object' must be of class 'master_matrix'")
    }
  }else {
    stop("Argument 'master_matrix' is necessary.")
  }

  cat("\n                     Summary of a master_matrix object\n")
  cat("---------------------------------------------------------------------------\n\n")
  cat("Data matrix summary:\n\n")
  print(summary(object$data_matrix))
  if (!is.null(object$preselected_sites)) {
    cat("\n\nSites preselected by user:\n\n")
    print(object$preselected_sites[, 1:3])
  }
  cat("\n\nRegion of interest:\n\n")
  print(object$region)
  if (!is.null(object$mask)) {
    cat("\n\nMask to area of interest:\n\n")
    print(object$mask)
  }
  if (!is.null(object$PCA_results)) {
    cat("\n\nPCA results:\n\n")
    print(object$PCA_results)
  }
}

