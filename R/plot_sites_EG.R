#' Representation of sites selected to be surveyed
#'
#' @description Plots representing sites (all and selected for survey) in
#' environmental and/or geographic space.
#'
#' @aliases plot_sites_EG plot_sites_E plot_sites_G
#'
#' @param master_selection master_selection object derived from functions
#' \code{\link{random_selection}}, \code{\link{uniformG_selection}},
#' \code{\link{uniformE_selection}}, or \code{\link{EG_selection}}.
#' @param selection_type (character) type of selection depending on the function
#' used to select sites. The options available are "random"
#' (\code{\link{random_selection}}), "G" (\code{\link{uniformG_selection}}),
#' "E" (\code{\link{uniformE_selection}}), and "EG"
#' (\code{\link{EG_selection}}).
#' @param variable_1 (character or numeric) name or position of the first
#' variable (x-axis) to be plotted in environmental space. Default = NULL,
#' required when \code{selection_type} = "random" or "G".
#' @param variable_2 (character or numeric) name or position of the second
#' variable (y-axis) to be plotted in environmental space. It must be different
#' from the first one. Default = NULL, required when \code{selection_type} =
#' "random" or "G".
#' @param selection_number (numeric) number of selection to be plotted.
#' Default = 1.
#' @param region_border (logical) whether to add region border to the plot.
#' Default = TRUE.
#' @param mask_border (logical) whether to add mask border to the plot. Ignored
#' if mask is not present in \code{master_selection}. Default = FALSE.
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
#' @param pch_all (numeric) integer specifying a symbol when plotting all
#' points. Default = 16.
#' @param pch_sites (numeric) integer specifying a symbol when plotting points
#' of selected sites. Default = 16.
#' @param pch_pre (numeric) integer specifying a symbol when plotting points
#' of preselected sites. Default = 16. Ignored if preselected sites are
#' not present in \code{master_selection}.
#' @param add_main (logical) whether or not to add fixed titles to the plot.
#' Default = TRUE. Titles added are "Environmental space" and "Geographic
#' space".
#'
#' @return
#' \code{plot_sites_EG} returns a two-panel plot showing the selected sites.
#' They are show in both spaces, geographic and environmental.
#'
#' \code{plot_sites_E} returns a plot of selected sites in environmental space.
#'
#' \code{plot_sites_G} returns a plot of selected sites in geographic space.
#'
#' @usage
#' plot_sites_EG(master_selection, selection_type, variable_1 = NULL,
#'               variable_2 = NULL, selection_number = 1,
#'               region_border = TRUE, mask_border = FALSE, col_all = NULL,
#'               col_sites = NULL, col_pre = NULL, cex_all = 0.7,
#'               cex_sites = 1, cex_pre = 1, pch_all = 16, pch_sites = 16,
#'               pch_pre = 16, add_main = TRUE)
#'
#' @export
#' @importFrom sp plot
#' @importFrom graphics layout par plot.new text title
#' @importFrom maps map
#'
#' @rdname plot_sites_EG
#'
#' @examples
#' # Data
#' data("m_selection", package = "biosurvey")
#'
#' # Plotting
#' plot_sites_EG(m_selection, selection_type = "E")
#' plot_sites_E(m_selection, selection_type = "E")
#' plot_sites_G(m_selection, selection_type = "E")


plot_sites_EG <- function(master_selection, selection_type, variable_1 = NULL,
                          variable_2 = NULL, selection_number = 1,
                          region_border = TRUE, mask_border = FALSE,
                          col_all = NULL, col_sites = NULL, col_pre = NULL,
                          cex_all = 0.7, cex_sites = 1, cex_pre = 1,
                          pch_all = 16, pch_sites = 16, pch_pre = 16,
                          add_main = TRUE) {
  # Initial tests
  if (missing(master_selection)) {
    stop("Argument 'master_selection' is required to produce the plot.")
  }
  if (class(master_selection)[1] != "master_selection") {
    stop("Object defined in 'master_selection' is not valid, see function's help.")
  }

  ## Par settings
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  # Plot
  ## Main-layout
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

  ## Environmental space
  plot_sites_E(master_selection, selection_type, variable_1, variable_2,
               selection_number, col_all, col_sites, col_pre, cex_all,
               cex_sites, cex_pre, pch_all, pch_sites, pch_pre)

  ## Geographic space
  par(mar = c(3.5, rep(0.5, 3)))
  plot_sites_G(master_selection, selection_type, selection_number,
               region_border, mask_border, col_all, col_sites, col_pre,
               cex_all, cex_sites, cex_pre, pch_all, pch_sites, pch_pre)
}



#' @rdname plot_sites_EG
#' @export
#' @param main (character) the main title for the plot.
#' @param xlab (character) label for the x axis. The default, NULL, uses
#' variable_1.
#' @param ylab (character) label for the y axis. The default, NULL, uses
#' variable_2.
#' @usage
#' plot_sites_E(master_selection, selection_type, variable_1 = NULL,
#'              variable_2 = NULL, selection_number = 1, col_all = NULL,
#'              col_sites = NULL, col_pre = NULL, cex_all = 0.7,
#'              cex_sites = 1, cex_pre = 1, pch_all = 16, pch_sites = 16,
#'              pch_pre = 16, main = "", xlab = NULL, ylab = NULL)

plot_sites_E <- function(master_selection, selection_type, variable_1 = NULL,
                         variable_2 = NULL, selection_number = 1,
                         col_all = NULL, col_sites = NULL, col_pre = NULL,
                         cex_all = 0.7, cex_sites = 1, cex_pre = 1,
                         pch_all = 16, pch_sites = 16, pch_pre = 16,
                         main = "", xlab = NULL, ylab = NULL) {
  # Initial tests
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

  # Preparing data
  evars <- c(variable_1, variable_2)

  # Colors
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

  # labels
  if (is.null(xlab)) {
    xlab <- variable_1
  }
  if (is.null(ylab)) {
    ylab <- variable_2
  }

  # Plot
  ## Environmental space
  plot(master_selection$data_matrix[, evars], col = col_all, pch = pch_all,
       cex = cex_all, bty = "l", xlab = "", ylab = "", main = main)
  title(xlab = xlab, line = 2.3, cex.lab = 1.1)
  title(ylab = ylab, line = 2.3, cex.lab = 1.1)

  ## Selected sites
  selected_data <- master_selection[[selection_type]][[selection_number]]
  points(selected_data[, evars], pch = pch_sites, cex = cex_sites,
         col = col_sites)

  ## Preselected sites
  precon <- sel_args$arguments$use_preselected_sites
  if (!is.null(master_selection$preselected_sites) & precon == TRUE) {
    points(master_selection$preselected_sites[, evars], pch = pch_pre,
           cex = cex_pre, col = col_pre)
  }

}


#' @rdname plot_sites_EG
#' @export
#' @usage
#' plot_sites_G(master_selection, selection_type, selection_number = 1,
#'              region_border = TRUE, mask_border = FALSE, col_all = NULL,
#'              col_sites = NULL, col_pre = NULL, cex_all = 0.7,
#'              cex_sites = 1, cex_pre = 1, pch_all = 16, pch_sites = 16,
#'              pch_pre = 16)

plot_sites_G <- function(master_selection, selection_type, selection_number = 1,
                         region_border = TRUE, mask_border = FALSE,
                         col_all = NULL, col_sites = NULL, col_pre = NULL,
                         cex_all = 0.7, cex_sites = 1, cex_pre = 1,
                         pch_all = 16, pch_sites = 16, pch_pre = 16) {
  # Initial tests
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
  }

  # Preparing data
  gvars <- c("Longitude", "Latitude")

  # Colors
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

  # Plot
  ## Box to plot
  boxpam <- t(master_selection$region@bbox)
  crsm <- master_selection$region@proj4string
  boxpam <- sp::SpatialPointsDataFrame(boxpam, data.frame(boxpam),
                                       proj4string = crsm)

  ## The plot
  sp::plot(boxpam, col = NA)
  maps::map(fill = TRUE, col = "gray97", lforce = "n",
            border = "gray80", add = TRUE)
  box(which = "plot")

  ## Region
  if (region_border == TRUE) {
    sp::plot(master_selection$region, border = "gray50", add = TRUE)
  }
  if (mask_border == TRUE & !is.null(master_selection$mask)) {
    sp::plot(master_selection$mask, border = "gray50", add = TRUE)
  }

  ## Selected sites
  selected_data <- master_selection[[selection_type]][[selection_number]]
  points(selected_data[, gvars], pch = pch_sites, cex = cex_sites,
         col = col_sites)

  ## Preselected sites
  precon <- sel_args$arguments$use_preselected_sites
  if (!is.null(master_selection$preselected_sites) & precon == TRUE) {
    points(master_selection$preselected_sites[, gvars], pch = pch_pre,
           cex = cex_pre, col = col_pre)
  }
}
