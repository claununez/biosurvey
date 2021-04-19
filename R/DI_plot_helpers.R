#' Plot dissimilarities among sets of selected sites as a dendrogram
#'
#' @param DI_selected_sites list of results obtained with function
#' \code{\link{selected_sites_DI}}.
#' @param labels (character) vector of labels for the tips of the tree.
#' The default, NULL, uses names of sets of selected sites. If labels = FALSE
#' no tip labels are plotted.
#' @param xlab (character) label for x axis of plot. Default = "".
#' @param ylab (character) label for y axis of plot. Default = "Distance".
#' @param main (character) title for the plot. Default = "Cluster dendrogram".
#' @param sub (character) subtitle for the plot. Plotted below the label of
#' the x axis.
#' @param ... other arguments to be passed to plot method for objects of class
#' "\code{hclust}". See more details in \code{\link[stats]{hclust}}.
#'
#' @return
#' A dendrogram plot of a "\code{hclust}" object.
#'
#' @usage
#' DI_dendrogram(DI_selected_sites, labels = NULL, xlab = "",
#'               ylab = "Distance", main = "Cluster dendrogram",
#'               sub = "", ...)
#'
#' @export
#' @import stats
#'
#' @examples
#' # Data
#' data("b_pam", package = "biosurvey")
#' data("m_selection", package = "biosurvey")
#'
#' # Subsetting base PAM according to selections
#' sub_pam_all <- subset_PAM(b_pam, m_selection, selection_type = "all")
#'
#' # Calculating dissimilarities
#' DI_sel <- selected_sites_DI(sub_pam_all)
#'
#' # Plot
#' DI_dendrogram(DI_sel)

DI_dendrogram <- function(DI_selected_sites, labels = NULL,
                          xlab = "", ylab = "Distance",
                          main = "Cluster dendrogram", sub = "", ...) {

  # Initial tests
  if (missing(DI_selected_sites)) {
    stop("Argument 'DI_selected_sites' must be defined.")
  }

  # Plot
  if (is.null(labels)) {
    labels <- gsub("_", " ", DI_selected_sites$cluster_selections$labels)
  }

  plot(DI_selected_sites$cluster_selections, labels = labels,
       ylab = ylab, xlab = xlab, main = main, las = 1, sub = sub, ...)
}
