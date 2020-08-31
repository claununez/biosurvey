#' Plot of PAM indices in geography
#'
#' @param PAM an object of class base_PAM.
#' @param index (character) code for the index to be plotted. Options are: "RI"
#' (Richness), "RIN" (Richness normalized), "DF" (Dispersion field), or "MCC"
#' (Mean composition covariance). Default = "RI".
#' @param master_selection a master_selection object derived from functions
#' \code{\link{random_selection}}, \code{\link{uniformG_selection}},
#' \code{\link{uniformE_selection}}, or \code{\link{EG_selection}}.
#' @param selection_type (character) Type of selection depending on the function
#' used to select sites. The options available are "random"
#' (\code{\link{random_selection}}), "G" (\code{\link{uniformG_selection}}),
#' "E" (\code{\link{uniformE_selection}}), and "EG" (\code{\link{EG_selection}}).
#' @param selection_number (numeric) number of selection to be plotted.
#' Default = 1.
#' @param col_pal color palette function to be used in defining colors for the
#' \code{index} to be plotted. The default, NULL, uses a color blind friendly
#' palette similar to viridis.
#' @param border color for cell borders of the PAM grid. The default, NULL, does
#' not plot any border.
#' @param col_sites color for selected sites. The default, NULL, uses
#' a red color to represent selected sites.
#' @param col_pre color for preselected sites. The default, NULL, uses
#' a purple color to represent preselected sites. Ignored if preselected sites are
#' not present in \code{master_selection}.
#' @param pch_sites (numeric) integer specifying a symbol when plotting points
#' of selected sites. Default = 16.
#' @param pch_pre (numeric) integer specifying a symbol when plotting points
#' of preselected sites. Default = 16. Ignored if preselected sites are
#' not present in \code{master_selection}.
#' @param cex (numeric) value by which plotting elements should be magnified
#' relative to the default. Default = 0.9
#'
#' @return
#' A plot of \code{index} represented in geography. Selected sites are added if
#' \code{master_selection} is defined.
#'
#' @usage
#' plot_PAM_geo(PAM, index = "RI", master_selection = NULL,
#'              selection_type = NULL, selection_number = 1,
#'              col_pal = NULL, border = NULL, col_sites = NULL,
#'              col_pre = NULL, pch_sites = 16, pch_pre = 16,cex = 0.9)
#'
#' @export
#' @importFrom sp plot
#' @importFrom graphics layout par plot.new
#' @importFrom maps map
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' # data
#' data("b_pam", package = "biosurvey")
#'
#' # plotting
#' plot_PAM_geo(b_pam, index = "RI")

plot_PAM_geo <- function(PAM, index = "RI", master_selection = NULL,
                         selection_type = NULL, selection_number = 1,
                         col_pal = NULL, border = NULL, col_sites = NULL,
                         col_pre = NULL, pch_sites = 16, pch_pre = 16,cex = 0.9) {
  if (missing(PAM)) {
    stop("Argument 'PAM' is missing.")
  }
  all_in <- c("RI", "RIN", "DF", "MCC")
  if (!index %in% all_in) {
    stop("Argument 'index' is not valid, options are: 'RI', 'RIN', 'DF', or 'MCC'.")
  }

  if (!is.null(master_selection)) {
    if (is.null(selection_type)) {
      stop("Argument 'selection_type' must be defined.")
    } else {
      if (!selection_type %in% c("random", "E", "G", "EG")) {
        stop("Argument 'selection_type' is not valid, options are: 'random', 'E', 'G', or 'EG'.")
      } else {
        selection_type <- paste0("selected_sites_", selection_type)
        sel_args <- attributes(master_selection[[selection_type]])
      }
    }

    # Where to visualize data
    where <- ifelse(!is.null(master_selection$mask), "mask", "region")
    gvars <- c("Longitude", "Latitude")
    precon <- sel_args$arguments$use_preselected_sites
  }

  # index selection
  PAM <- PAM_indices(PAM, indices = "all")
  g_indices <- names(PAM$PAM_indices)[-c(1, 3, 5, 7, 9:11)]
  names(g_indices) <- all_in

  # color definition
  if (is.null(col_pal)) {
    col_pal <- colorRampPalette(rev(c("#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb",
                                      "#41b6c4", "#1d91c0", "#225ea8", "#253494",
                                      "#081d58")))
  }
  if (is.null(border)) {
    border <- NA
  }

  if (is.null(col_sites) & is.null(col_pre)) {
    col_sites <- "#EC2F06"
    col_pre <- "#DD00FF"
  } else {
    if (is.null(col_sites)) {
      col_sites <- "#871B04"
    }
    if (is.null(col_pre)) {
      col_pre <- "#DD00FF"
    }
  }

  rfactor <- range(PAM$PAM_indices[[g_indices[index]]])
  ifactor <- as.factor(PAM$PAM_indices[[g_indices[index]]])
  n <- length(levels(ifactor))
  col <- col_pal(n)

  # par settings
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  # plotting
  layout(matrix(1:2, 1, byrow = T), widths = c(10, 1.5))
  par(cex = cex, mar = rep(0, 4))

  sp::plot(PAM$PAM, border = "transparent")
  maps::map(fill = TRUE, col = "gray97", lforce = "n",
            border = "gray80", add = TRUE)
  sp::plot(PAM$PAM, col = col[ifactor], border = border, add = TRUE)
  box()

  if (!is.null(master_selection)) {
    if (is.null(master_selection$mask)) {
      sp::plot(master_selection$region, border = "gray70", add = TRUE)
    }
    sp::plot(master_selection[[where]], border = "gray60", add = TRUE)

    ## selected sites
    selected_data <- master_selection[[selection_type]][[selection_number]]
    points(selected_data[, gvars], pch = pch_sites, col = col_sites)

    ## preselected sites
    if (!is.null(master_selection$preselected_sites) & precon == TRUE) {
      points(master_selection$preselected_sites[, gvars], pch = pch_pre,
             col = col_pre)
    }
  }

  plot.new()
  bar_legend(rfactor, col = col, title = gsub("_", " ", g_indices[index]),
             round = 3, label_x = 0.5, labels_y = c(0.18, 0.87))
}



