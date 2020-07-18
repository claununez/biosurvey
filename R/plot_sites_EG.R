#' Representation of sites selected to be surveyed
#'
#' @description Creates a two-panel plot representing sites (all and selected
#' for survey) in both spaces, environmental and geographic.
#'
#' @param master_selection a master_selection object derived from functions
#' \code{\link{random_selection}}, \code{\link{uniformG_selection}},
#' \code{\link{uniformE_selection}}, or \code{\link{EG_selection}}.
#' @param variable_1 (character or numeric) name or position of the first
#' variable (X axis).
#' @param variable_2 (character or numeric) name or position of the second
#' variable (Y axis). It must be different from the first one.
#' @param selection_type (character) Type of selection depending on the function
#' used to select sites. The options available are "selected_sites_random"
#' (\code{\link{random_selection}}), "selected_sites_G"
#' (\code{\link{uniformG_selection}}), "selected_sites_E"
#' (\code{\link{uniformE_selection}}), and "selected_sites_EG"
#' (\code{\link{EG_selection}}).
#' @param selection_number (numeric) number of selection to be plotted.
#' Default = 1.
#' @param col_all colors for points in all points in the region of interest.
#' The default, NULL, uses a light gray color.
#' @param col_sites color for selected sites. The default, NULL, uses
#' a blue color to represent selected sites.
#' @param cex_all (numeric) value defining magnification of all points
#' relative to the default. Default = 0.7.
#' @param cex_sites (numeric) value defining magnification of selected sites
#' relative to the default. Default = 1.
#' @param pch_all (numeric) integer specifying a symbol when plotting all points.
#' Default = 16.
#' @param pch_sites (numeric) integer specifying a symbol when plotting points
#' of selected sites. Default = 16.
#' @param add_main (logical) whether or not to add fixed titles to the plot.
#' Default = TRUE. Titles added are "Environmental space" and "Geographic space".
#'
#' @return
#' A two-panel plot showing the selected sites. They are show in both spaces,
#' geographic and environmental.
#'
#' @usage
#' plot_sites_EG(master_selection, variable_1, variable_2, selection_type,
#'               selection_number = 1, col_all = NULL, col_sites = NULL,
#'               cex_all = 0.7, cex_sites = 1, pch_all = 16,
#'               pch_sites = 16, add_main = TRUE)
#'
#' @export
#' @importFrom sp plot
#' @importFrom graphics layout par plot.new text title
#' @importFrom maps map
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#'
#' # Making blocks for analysis
#' m_blocks <- make_blocks(m_matrix, variable_1 = "PC1",
#'                         variable_2 = "PC2", n_cols = 10, n_rows = 10,
#'                         block_type = "equal_area")
#'
#' # Checking column names
#' colnames(m_blocks$master_matrix)
#'
#' # Selecting sites uniformly in E space
#' selectionE <- uniformE_selection(m_blocks, variable_1 = "PC1", variable_2 = "PC2",
#'                                  selection_from = "block_centroids",
#'                                  expected_points = 25, max_n_samplings = 1,
#'                                  initial_distance = 1, increase = 0.1,
#'                                  replicates = 5, set_seed = 1)
#'
#' # Plotting
#' plot_sites_EG(selectionE, variable_1 = "PC1", variable_2 = "PC2",
#'               selection_type = "selected_sites_E")


plot_sites_EG <- function(master_selection, variable_1, variable_2, selection_type,
                          selection_number = 1, col_all = NULL, col_sites = NULL,
                          cex_all = 0.7, cex_sites = 1, pch_all = 16,
                          pch_sites = 16, add_main = TRUE) {
  # initial tests
  if (missing(master_selection)) {
    stop("Argument 'master_selection' is required to produce the plot.")
  }
  if (missing(variable_1)) {
    stop("Argument 'variable_1' is required to produce the plot.")
  }
  if (missing(variable_2)) {
    stop("Argument 'variable_2' is required to produce the plot.")
  }
  if (class(master_selection)[1] != "master_selection") {
    stop("Object defined in 'master_selection' is not valid, see function's help.")
  }
  if (missing(selection_type)) {
    stop("Argument 'selection_type' is required to produce the plot.")
  }
  if (!selection_type %in% c("selected_sites_random", "selected_sites_E",
                             "selected_sites_G", "selected_sites_EG")) {
    stop("Argument 'selection_type' is not valid, options are:\n'selected_sites_random', 'selected_sites_E', 'selected_sites_G', or 'selected_sites_EG'.")
  }

  # Where to visualize data
  where <- ifelse(!is.null(master_selection$mask), "mask", "region")

  # preparing data
  xlim <- master_selection[[where]]@bbox[1, ]
  ylim <- master_selection[[where]]@bbox[2, ]

  gvars <- c("Longitude", "Latitude")
  evars <- c(variable_1, variable_2)

  # colors
  if (is.null(col_all) & is.null(col_sites)) {
    col_all <- "#E1E1E1"
    col_sites <- "#3B22CB"
  } else {
    if (is.null(col_all)) {
      col_all <- "#E1E1E1"
    }
    if (is.null(col_sites)) {
      col_sites <- "#3B22CB"
    }
  }

  ## par settings
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  # plot
  ## main-layout
  if (add_main == TRUE) {
    layout(matrix(1:4, 2, byrow = T), widths = c(10, 10), heights = c(1, 10))
    par(mar = rep(0, 4), cex = 0.7)

    plot.new()
    text(0.5, 0.5, "Environmental space", cex = 1.3)
    plot.new()
    text(0.5, 0.5, "Geographic space", cex = 1.3)

    par(mar = c(3.5, 3.5, 0.5, 1))
  } else {
    if (all.equal(par()$mfrow, c(1, 1)) == TRUE) {
      par(mfrow = c(1, 2))
    }
    par(mar = c(3.5, 3.5, 0.5, 1), cex = 0.7)
  }

  ## environmental space
  plot(master_selection$master_matrix[, evars], col = col_all, pch = pch_all,
       cex = cex_all, bty = "l", xlab = "", ylab = "")
  title(xlab = variable_1, line = 2.3, cex.lab = 1.1)
  title(ylab = variable_2, line = 2.3, cex.lab = 1.1)

  ## selected sites
  selected_data <- master_selection[[selection_type]][[selection_number]]
  points(selected_data[, evars], pch = pch_sites, cex = cex_sites, col = col_sites)

  ## geographic space
  par(mar = rep(0.5, 4))
  sp::plot(master_selection[[where]], border = "transparent")
  maps::map(fill = TRUE, col = "gray97", lforce = "n",
            border = "gray80", add = TRUE)
  box(which = "plot")
  points(master_selection$master_matrix[, gvars], pch = pch_all, cex = cex_all,
         col = col_all)
  if (is.null(master_selection$mask)) {
    sp::plot(master_selection$region, border = "gray70", add = TRUE)
  }
  sp::plot(master_selection[[where]], add = TRUE)

  ## selected sites
  points(selected_data[, gvars], pch = pch_sites, cex = cex_sites, col = col_sites)
}
