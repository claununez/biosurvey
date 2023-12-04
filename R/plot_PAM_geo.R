#' Plot of PAM indices in geography
#'
#' @param PAM object of class base_PAM.
#' @param index (character) code for the index to be plotted. Options are: "RI"
#' (Richness), "RIN" (Richness normalized), "DF" (Dispersion field), or "MCC"
#' (Mean composition covariance). Default = "RI".
#' @param master_selection master_selection object derived from functions
#' \code{\link{random_selection}}, \code{\link{uniformG_selection}},
#' \code{\link{uniformE_selection}}, or \code{\link{EG_selection}}.
#' @param region_border (logical) whether to add region border to the plot.
#' Default = TRUE.
#' @param mask_border (logical) whether to add mask border to the plot. Ignored
#' if mask is not present in \code{master_selection}. Default = FALSE.
#' @param selection_type (character) Type of selection depending on the function
#' used to select sites. The options available are "random"
#' (\code{\link{random_selection}}), "G" (\code{\link{uniformG_selection}}),
#' "E" (\code{\link{uniformE_selection}}), and "EG" (\code{\link{EG_selection}}).
#' @param selection_number (numeric) number of selection to be plotted.
#' Default = 1.
#' @param col_pal color palette function to be used in defining colors for the
#' \code{index} to be plotted. The default, NULL, uses a simple color blind
#' friendly palette similar to viridis, \code{\link{purplow}}. Warmer colors
#' indicate higher values.
#' @param border color for cell borders of the PAM grid. The default, NULL, does
#' not plot any border.
#' @param col_sites color for selected sites. The default, NULL, uses
#' a red color to represent selected sites.
#' @param col_pre color for preselected sites. The default, NULL, uses
#' a purple color to represent preselected sites. Ignored if preselected sites
#' are not present in \code{master_selection}.
#' @param pch_sites (numeric) integer specifying a symbol when plotting points
#' of selected sites. Default = 16.
#' @param pch_pre (numeric) integer specifying a symbol when plotting points
#' of preselected sites. Default = 16. Ignored if preselected sites are
#' not present in \code{master_selection}.
#'
#' @return
#' A plot of \code{index} represented in geography. Selected sites are added if
#' \code{master_selection} is defined.
#'
#' @usage
#' plot_PAM_geo(PAM, index = "RI", master_selection = NULL,
#'              region_border = TRUE, mask_border = FALSE,
#'              selection_type = NULL, selection_number = 1,
#'              col_pal = NULL, border = NULL, col_sites = NULL,
#'              col_pre = NULL, pch_sites = 16, pch_pre = 16)
#'
#' @export
#' @importFrom sp plot
#' @importFrom graphics layout par plot.new
#' @importFrom maps map
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' # Data
#' b_pam <- read_PAM(system.file("extdata/b_pam.rds",
#'                               package = "biosurvey"))
#'
#' # Plotting
#' plot_PAM_geo(b_pam, index = "RI")
#'
#' # You can add a legend with
#' legend_bar(position = "bottomleft", col = purplow(8), title = "Richness",
#'            labels = c("Low", "High"))


plot_PAM_geo <- function(PAM, index = "RI", master_selection = NULL,
                         region_border = TRUE, mask_border = FALSE,
                         selection_type = NULL, selection_number = 1,
                         col_pal = NULL, border = NULL, col_sites = NULL,
                         col_pre = NULL, pch_sites = 16, pch_pre = 16) {
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
    gvars <- c("Longitude", "Latitude")
    precon <- sel_args$arguments$use_preselected_sites
  }

  # Index selection
  inPAM_in <- names(PAM$PAM_indices)

  if (!index %in% inPAM_in) {
    PAM <- PAM_indices(PAM, indices = "all")
  }

  g_indices <- names(PAM$PAM_indices)[-c(1, 3, 5, 7, 9:11)]
  names(g_indices) <- all_in

  # Color definition
  if (is.null(col_pal)) {
    col_pal <- purplow
  }
  border <- ifelse(is.null(border), NA, border)

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

  ifactor <- as.factor(PAM$PAM_indices[[g_indices[index]]])
  n <- length(levels(ifactor))
  col <- col_pal(n)

  # Box to plot
  boxpam <- matrix(terra::ext(PAM$PAM), nrow = 2)
  boxpam <- terra::vect(boxpam, crs = terra::crs(PAM$PAM))

  # Plotting
  terra::plot(boxpam, col = NA, axes = FALSE, legend = FALSE)
  maps::map(fill = TRUE, col = "gray97", lforce = "n",
            border = "gray80", add = TRUE)
  terra::plot(PAM$PAM, col = col[ifactor], border = border, add = TRUE)
  box()

  if (!is.null(master_selection)) {
    if (region_border == TRUE) {
      terra::plot(master_selection$region, border = "gray50", add = TRUE)
    }
    if (mask_border == TRUE & !is.null(master_selection$mask)) {
      terra::plot(master_selection$mask, border = "gray50", add = TRUE)
    }

    ## Selected sites
    selected_data <- master_selection[[selection_type]][[selection_number]]
    points(selected_data[, gvars], pch = pch_sites, col = col_sites)

    ## Preselected sites
    if (!is.null(master_selection$preselected_sites) & precon == TRUE) {
      points(master_selection$preselected_sites[, gvars], pch = pch_pre,
             col = col_pre)
    }
  }
}
