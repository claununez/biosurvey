#' Plotting dissimilarity indices withing and among sets of selected sites
#'
#' @description Creates matrix-like plots of dissimilarities found among
#' communities of species in distinct sites selected or sets of sites selected
#' for sampling.
#'
#' @param DI_selected_sites list of results obtained with function
#' \code{\link{selected_sites_DI}}.
#' @param selection_type type of selection to be considered when creating DI
#' matrix plot. Options are: "selections", "random", "E", "G", and "EG".
#' The default, "selections", plots a comparison among all selection types.
#' @param selection_number (numeric) number of selection to be plotted.
#' Default = 1. Ignored if \code{selection_type} = "selections".
#' @param values (logical) whether or not to add values of dissimilarity.
#' Default = TRUE.
#' @param col a list of colors derived from a palette. Default =
#' \code{heat.colors(12, rev = TRUE)}.
#' @param xlab (character) label for x axis of plot. Default = "Number of sites".
#' @param ylab (character) label for y axis of plot. Default = "Species".
#'
#' @return
#' A plot of a matrix of dissimilarities among sites selected for sampling, or
#' among sets of sampling sites selected. Random is abbreviated as "R" in labels.
#'
#' @usage
#' plot_DI(DI_selected_sites, selection_type = "selections",
#'         selection_number = 1, values = TRUE,
#'         col = heat.colors(12, rev = TRUE),
#'         xlab = "", ylab = "")
#'
#' @export
#' @importFrom graphics image axis text box
#' @importFrom grDevices heat.colors
#'
#' @examples
#' # Data
#' b_pam <- read_PAM(system.file("extdata/b_pam.rds",
#'                               package = "biosurvey"))
#' m_selection <- read_master(system.file("extdata/m_selection.rds",
#'                                        package = "biosurvey"))
#'
#' # Subsetting base PAM according to selections
#' sub_pam_all <- subset_PAM(b_pam, m_selection, selection_type = "all")
#'
#' # Calculating dissimilarities
#' DI_sel <- selected_sites_DI(sub_pam_all)
#'
#' # Plotting
#' plot_DI(DI_sel)

plot_DI <- function(DI_selected_sites, selection_type = "selections",
                    selection_number = 1, values = TRUE,
                    col = heat.colors(12, rev = TRUE),
                    xlab = "", ylab = "") {

  # Initial tests
  if (missing(DI_selected_sites)) {
    stop("Argument 'DI_selected_sites' must be defined.")
  }
  if (!selection_type %in% c("selections", "random", "E", "G", "EG")) {
    stop("Argument 'selection_type' is not valid, options are: 'selections', 'random', 'E', 'G', or 'EG'.")
  }

 # Preparing data for plotting
  if (selection_type == "selections") {
    selection_type <- paste0("DI_", selection_type)
    mat <- data.matrix(DI_selected_sites[[selection_type]])
  } else {
    selection_type <- paste0("DI_selected_sites_", selection_type)
    mat <- data.matrix(DI_selected_sites[[selection_type]][[selection_number]])
  }

  cnam <- colnames(mat)
  cnam <- gsub("_", " ", gsub("random", "R", cnam))
  mat <- apply(mat, 2, rev)
  cmat <- sprintf("%0.1f", t(mat))
  dm <- ncol(mat)

  # Plotting
  image(1:dm, 1:dm, t(mat), col = col, axes = FALSE, xlab = "", ylab = "")

  axis(3, 1:dm, cnam, cex.axis = 0.6, las = 3)
  axis(2, 1:dm, rev(cnam), cex.axis = 0.6, las = 1)
  box()

  if (values == TRUE) {
    text(expand.grid(1:dm, 1:dm), cmat, cex = 0.7)
  }
}





#' Plot dissimilarities withing and among sets of selected sites as a dendrogram
#'
#' @param DI_selected_sites list of results obtained with function
#' \code{\link{selected_sites_DI}}.
#' @param selection_type type of selection to be considered when creating DI
#' matrix plot. Options are: "selections", "random", "E", "G", and "EG".
#' The default, "selections", plots a comparison among all selection types.
#' @param selection_number (numeric) number of selection to be plotted.
#' Default = 1. Ignored if \code{selection_type} = "selections".
#' @param labels (character) vector of labels for the tips of the tree.
#' The default, NULL, uses names of sets of selected sites. If labels = FALSE
#' no tip labels are plotted.
#' @param xlab (character) label for x-axis of plot. Default = "".
#' @param ylab (character) label for y-axis of plot. Default = "Distance".
#' @param main (character) title for the plot. Default = "Cluster dendrogram".
#' @param sub (character) subtitle for the plot. Plotted below the label of
#' the x-axis.
#' @param ... other arguments to be passed to plot method for objects of class
#' "\code{hclust}". See more details in \code{\link[stats]{hclust}}.
#'
#' @return
#' A dendrogram plot of a "\code{hclust}" object.
#'
#' @usage
#' DI_dendrogram(DI_selected_sites, selection_type = "selections",
#'               selection_number = 1, labels = NULL, xlab = "",
#'               ylab = "Distance", main = "Cluster dendrogram",
#'               sub = "", ...)
#'
#' @export
#' @import stats
#'
#' @examples
#' # Data
#' b_pam <- read_PAM(system.file("extdata/b_pam.rds",
#'                               package = "biosurvey"))
#' m_selection <- read_master(system.file("extdata/m_selection.rds",
#'                                        package = "biosurvey"))
#'
#' # Subsetting base PAM according to selections
#' sub_pam_all <- subset_PAM(b_pam, m_selection, selection_type = "all")
#'
#' # Calculating dissimilarities
#' DI_sel <- selected_sites_DI(sub_pam_all)
#'
#' # Plot
#' DI_dendrogram(DI_sel)

DI_dendrogram <- function(DI_selected_sites, selection_type = "selections",
                          selection_number = 1, labels = NULL,
                          xlab = "", ylab = "Distance",
                          main = "Cluster dendrogram", sub = "", ...) {

  # Initial tests
  if (missing(DI_selected_sites)) {
    stop("Argument 'DI_selected_sites' must be defined.")
  }
  if (!selection_type %in% c("selections", "random", "E", "G", "EG")) {
    stop("Argument 'selection_type' is not valid, options are: 'selections'', 'random', 'E', 'G', or 'EG'.")
  }

  # Preparing data for plotting
  stype <- paste0("cluster_", selection_type)
  if (selection_type == "selections") {
    csel <- DI_selected_sites[[stype]]
  } else {
    csel <- DI_selected_sites[[stype]][[selection_number]]
  }

  # Plot
  if (is.null(labels)) {
    labels <- gsub("_", " ", csel$labels)
  }

  plot(csel, labels = labels, ylab = ylab, xlab = xlab, main = main,
       sub = sub, ...)
}
