#' Representation of sites selected to be surveyed
#'
#' @description Creates a two-panel plot representing sites (all and selected
#' for survey) in both spaces, environmental and geographic.
#'
#' @param master_selection a master_selection object derived from functions
#' \code{\link{random_selection}}, \code{\link{uniformG_selection}},
#' \code{\link{uniformE_selection}}, or \code{\link{EG_selection}}.
#' @param selection_type (character) Type of selection depending on the function
#' used to select sites. The options available are "random"
#' (\code{\link{random_selection}}), "G" (\code{\link{uniformG_selection}}),
#' "E" (\code{\link{uniformE_selection}}), and "EG" (\code{\link{EG_selection}}).
#' @param variable_1 (character or numeric) name or position of the first
#' variable (X axis) to be plotted in environmental space. Default = NULL,
#' required when \code{selection_type} = "random" or "G".
#' @param variable_2 (character or numeric) name or position of the second
#' variable (Y axis) to be plotted in environmental space. It must be different
#' from the first one. Default = NULL, required when \code{selection_type} =
#' "random" or "G".
#' @param selection_number (numeric) number of selection to be plotted.
#' Default = 1.
#' @param col_all colors for points in all points in the region of interest.
#' The default, NULL, uses a light gray color.
#' @param col_sites color for selected sites. The default, NULL, uses
#' a blue color to represent selected sites.
#' @param col_pre color for preselected sites. The default, NULL, uses
#' a red color to represent preselected sites. Ignored if preselected sites are
#' not present in \code{master_selection}.
#' @param cex_all (numeric) value defining magnification of all points
#' relative to the default. Default = 0.7.
#' @param cex_sites (numeric) value defining magnification of selected sites
#' relative to the default. Default = 1.
#' @param cex_pre (numeric) value defining magnification of preselected sites
#' relative to the default. Default = 1. Ignored if preselected sites are
#' not present in \code{master_selection}.
#' @param pch_all (numeric) integer specifying a symbol when plotting all points.
#' Default = 16.
#' @param pch_sites (numeric) integer specifying a symbol when plotting points
#' of selected sites. Default = 16.
#' @param pch_pre (numeric) integer specifying a symbol when plotting points
#' of preselected sites. Default = 16. Ignored if preselected sites are
#' not present in \code{master_selection}.
#' @param add_main (logical) whether or not to add fixed titles to the plot.
#' Default = TRUE. Titles added are "Environmental space" and "Geographic space".
#'
#' @return
#' A two-panel plot showing the selected sites. They are show in both spaces,
#' geographic and environmental.
#'
#' @usage
#' plot_sites_EG(master_selection, selection_type, variable_1 = NULL,
#'               variable_2 = NULL, selection_number = 1, col_all = NULL,
#'               col_sites = NULL, col_pre = NULL, cex_all = 0.7,
#'               cex_sites = 1, cex_pre = 1, pch_all = 16, pch_sites = 16,
#'               pch_pre = 16, add_main = TRUE)
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
#' colnames(m_blocks$data_matrix)
#'
#' # Selecting sites uniformly in E space
#' selectionE <- uniformE_selection(m_blocks, variable_1 = "PC1", variable_2 = "PC2",
#'                                  selection_from = "block_centroids",
#'                                  expected_points = 25, max_n_samplings = 1,
#'                                  initial_distance = 1, increase = 0.1,
#'                                  replicates = 5, set_seed = 1)
#'
#' # Plotting
#' plot_sites_EG(selectionE, selection_type = "E")


plot_sites_EG <- function(master_selection, selection_type, variable_1 = NULL,
                          variable_2 = NULL, selection_number = 1, col_all = NULL,
                          col_sites = NULL, col_pre = NULL, cex_all = 0.7,
                          cex_sites = 1, cex_pre = 1, pch_all = 16, pch_sites = 16,
                          pch_pre = 16, add_main = TRUE) {
  # initial tests
  if (missing(master_selection)) {
    stop("Argument 'master_selection' is required to produce the plot.")
  }
  if (class(master_selection)[1] != "master_selection") {
    stop("Object defined in 'master_selection' is not valid, see function's help.")
  }
  if (missing(selection_type)) {
    stop("Argument 'selection_type' is required to produce the plot.")
  }

  if (!selection_type %in% c("random", "E", "G", "EG")) {
    stop("Argument 'selection_type' is not valid, options are: 'random', 'E', 'G', or 'EG'.")
  } else {
    selection_type <- paste0("selected_sites_", selection_type)
    sel_args <- attributes(master_selection[[selection_type]])

    if (is.null(master_selection[[selection_type]])) {
      stop("'selection_type' defined is NULL in 'master_selection'.")
    }

    if (selection_type %in% c("selected_sites_random", "selected_sites_G")) {
      if (missing(variable_1)) {
        stop("Argument 'variable_1' is required to produce the plot.")
      }
      if (missing(variable_2)) {
        stop("Argument 'variable_2' is required to produce the plot.")
      }
    } else {
      variable_1 <- sel_args$arguments$variable_1
      variable_2 <- sel_args$arguments$variable_2
    }
  }

  # Where to visualize data
  where <- ifelse(!is.null(master_selection$mask), "mask", "region")

  # preparing data
  xlim <- master_selection[[where]]@bbox[1, ]
  ylim <- master_selection[[where]]@bbox[2, ]

  gvars <- c("Longitude", "Latitude")
  evars <- c(variable_1, variable_2)

  # colors
  if (is.null(col_all) & is.null(col_sites) & is.null(col_pre)) {
    col_all <- "#E1E1E1"
    col_sites <- "#3B22CB"
    col_pre <- "#EC2F06"
  } else {
    if (is.null(col_all)) {
      col_all <- "#E1E1E1"
    }
    if (is.null(col_sites)) {
      col_sites <- "#3B22CB"
    }
    if (is.null(col_pre)) {
      col_pre <- "#EC2F06"
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
  plot(master_selection$data_matrix[, evars], col = col_all, pch = pch_all,
       cex = cex_all, bty = "l", xlab = "", ylab = "")
  title(xlab = variable_1, line = 2.3, cex.lab = 1.1)
  title(ylab = variable_2, line = 2.3, cex.lab = 1.1)

  ## selected sites
  selected_data <- master_selection[[selection_type]][[selection_number]]
  points(selected_data[, evars], pch = pch_sites, cex = cex_sites, col = col_sites)

  ## preselected sites
  precon <- sel_args$arguments$use_preselected_sites
  if (!is.null(master_selection$preselected_sites) & precon == TRUE) {
    points(master_selection$preselected_sites[, evars], pch = pch_pre,
           cex = cex_pre, col = col_pre)
  }

  ## geographic space
  par(mar = rep(0.5, 4))
  sp::plot(master_selection[[where]], border = "transparent")
  maps::map(fill = TRUE, col = "gray97", lforce = "n",
            border = "gray80", add = TRUE)
  box(which = "plot")

  ## region
  if (is.null(master_selection$mask)) {
    sp::plot(master_selection$region, border = "gray70", add = TRUE)
  }
  sp::plot(master_selection[[where]], col = col_all, add = TRUE)

  ## selected sites
  points(selected_data[, gvars], pch = pch_sites, cex = cex_sites, col = col_sites)

  ## preselected sites
  if (!is.null(master_selection$preselected_sites) & precon == TRUE) {
    points(master_selection$preselected_sites[, gvars], pch = pch_pre,
           cex = cex_pre, col = col_pre)
  }
}
