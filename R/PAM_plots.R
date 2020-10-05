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
#' a purple color to represent preselected sites. Ignored if preselected sites
#' are not present in \code{master_selection}.
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





#' Plot new diversity-range diagram
#'
#' @param PAM_CS an object of class PAM_CS or a base_PAM object containing
#' a PAM_CS object as part of PAM_indices. These objects can be obtained using
#' the function \code{\link{prepare_PAM_CS}}.
#' @param add_significant (logical) whether to add statistically significant
#' values using a different symbol. Default = FALSE. If TRUE and values
#' indicating significance are not in \code{PAM_CS}, a message will be printed.
#' @param add_random_values (logical) whether to add values resulted from
#' the randomization process done when preparing \code{PAM_CS}. Default = FALSE.
#' Valid only if \code{add_significant} = TRUE, and randomized values are
#' present in \code{PAM_CS}.
#' @param col_all color code or name for all values. Default = "#8C8C8C".
#' @param col_significant_low color code or name for significant values below
#' confidence limits of random expectations. Default = "#000000".
#' @param col_significant_high color code or name for significant values above
#' confidence limits of random expectations. Default = "#000000".
#' @param col_random_values color code or name for randomized values.
#' Default = "D2D2D2".
#' @param pch_all point symbol to be used for all values. Default = 1.
#' @param pch_significant_low point symbol to be used for significant values
#' below confidence limits of random expectations. Default = 19.
#' @param pch_significant_high point symbol to be used for significant values
#' above confidence limits of random expectations. Default = 19.
#' @param pch_random_values point symbol to be used for randomized values.
#' Default = 1.
#' @param main main title for the plot. Default = NULL.
#' @param xlab label for the x axis. Default = NULL.
#' @param ylab label for the y axis. Default = NULL.
#' @param xlim x limits of the plot (x1, x2). The default, NULL, uses the range
#' of normalized richness.
#' @param ylim y limits of the plot. The default, NULL, uses the range of the
#' normalized values of the dispersion field. The second limit is increased
#' by adding the result of multiplying it by \code{ylim_expansion}, if
#' \code{add_legend} = TRUE.
#' @param ylim_expansion value used or expanding the \code{ylim}. Default = 0.25.
#' @param add_legend (logical) whether to add a legend describing information
#' relevant for interpreting the diagram. Default = TRUE.
#'
#' @return
#' A diversity-range plot with values of normalized richness in the x axis,
#' and normalized values of the dispersion field index divided by number of
#' species in the y axis.
#'
#' @usage
#' plot_PAM_CS(PAM_CS, add_significant = FALSE,
#'             add_random_values = FALSE, col_all = "#8C8C8C",
#'             col_significant_low = "#000000",
#'             col_significant_high = "#000000",
#'             col_random_values = "#D2D2D2", pch_all = 1,
#'             pch_significant_low = 19, pch_significant_high = 19,
#'             pch_random_values = 1, main = NULL,
#'             xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
#'             ylim_expansion = 0.25, add_legend = TRUE)
#'
#' @export
#' @importFrom graphics legend polygon
#'
#' @examples
#' # data
#' data("b_pam", package = "biosurvey")
#'
#' # preparing data for CS diagram
#' pcs <- prepare_PAM_CS(PAM = b_pam)
#'
#' # plot
#' plot_PAM_CS(pcs)

plot_PAM_CS <- function(PAM_CS, add_significant = FALSE,
                        add_random_values = FALSE, col_all = "#8C8C8C",
                        col_significant_low = "#000000",
                        col_significant_high = "#000000",
                        col_random_values = "#D2D2D2", pch_all = 1,
                        pch_significant_low = 19, pch_significant_high = 19,
                        pch_random_values = 1, main = NULL,
                        xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
                        ylim_expansion = 0.25, add_legend = TRUE) {

  if (!class(PAM_CS)[1] %in% c("base_PAM", "PAM_CS")) {
    stop("Class of 'PAM_CS' is not supported, see function's help.")
  }

  # preparing data
  if (class(PAM_CS)[1] == "base_PAM") {
    bp <- PAM_CS$PAM
    PAM <- PAM_CS$PAM_indices$CS_diagram
  } else {
    PAM <- PAM_CS
  }

  # Preparing values to be plotted or added
  s <- PAM$Species
  alfas <- PAM$Richness_normalized
  n <- PAM$Sites_cells
  fists <- PAM$Dispersion_field_normalized / s
  betty <- PAM$Beta_W

  # prepare vertex of plot limits
  vx <- PAM$Theoretical_boundaries$x
  vy <- PAM$Theoretical_boundaries$y
  sper <- round(PAM$Spearman_cor, 3)

  # plot elements
  if (is.null(main)) {
    main <- "Diversity-range plot"
  }
  if (is.null(xlab)) {
    xlab <- "Normalized richness"
  }
  if (is.null(ylab)) {
    ylab <- "Normalized dispersion field / S"
  }
  if (is.null(xlim)) {
    xlim <- range(vx)
  }
  if (is.null(ylim)) {
    add <- ifelse(add_legend == TRUE, ylim_expansion, 0)
    ylim <- range(vy) + c(0, range(vy)[2] * add)
  }

  # plot
  plot(alfas, fists, col = col_all, pch = pch_all, xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab, main = main)
  polygon(vx, vy, border = "#474747")
  if (add_legend == TRUE) {
    legend("topleft", bty = "n", inset = -0.02,
           legend = c(paste("N species (S) =", s), paste("N sites-cells =", n),
                      paste("Beta W =", round(betty, 3)),
                      as.expression(bquote("Spearman's" ~ r[s] ~ "=" ~ .(sper)))))
  }

  # randomized values
  if (add_random_values == TRUE) {
    if (!all(is.na(PAM$Randomized_DF))) {
      for (i in 1:ncol(PAM$Randomized_DF)) {
        points(alfas, PAM$Randomized_DF[, i], col = col_random_values,
               pch = pch_random_values, cex = 0.8)
      }
      points(alfas, fists, col = col_all, pch = pch_all)
    } else {
      message("Values from randomization process are missing in 'PAM_CS'")
    }
  }

  # significant values
  if (add_significant == TRUE) {
    if (!all(is.na(PAM$S_significance_id))) {
      sig_vals <- cbind(alfas, fists)[PAM$S_significance_id == 1, ]
      points(sig_vals, pch = pch_significant_low, col = col_significant_low)

      sig_vals <- cbind(alfas, fists)[PAM$S_significance_id == 2, ]
      points(sig_vals, pch = pch_significant_high, col = col_significant_high)
    } else {
      message("Values that indicate significance are missing in 'PAM_CS'")
    }
  }
}
