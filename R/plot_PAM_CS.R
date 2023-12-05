#' Representations of diversity and dispersion indices
#'
#' @description Graphic representations of results from
#' \code{\link{prepare_PAM_CS}}. Plots present the new range-diversity diagram
#' and geographic views of results. Geographic representations are only possible
#' when significant analyses were performed.
#'
#' @aliases plot_PAM_CS plot_PAM_CS_geo
#'
#' @param PAM_CS object of class PAM_CS or a base_PAM object containing
#' a PAM_CS object as part of PAM_indices. These objects can be obtained using
#' the function \code{\link{prepare_PAM_CS}}.
#' @param add_significant (logical) whether to add statistically significant
#' values using a different symbol. Default = FALSE. If TRUE and values
#' indicating significance are not in \code{PAM_CS}, a message will be printed.
#' @param add_random_values (logical) whether to add values resulted from
#' the randomization process done when preparing \code{PAM_CS}. Default = FALSE.
#' Valid only if \code{add_significant} = TRUE, and randomized values are
#' present in \code{PAM_CS}.
#' @param col_all color code or name for all values or those that are not
#' statistically significant. Default = "#CACACA".
#' @param col_significant_low color code or name for significant values below
#' confidence limits of random expectations. Default = "#000000".
#' @param col_significant_high color code or name for significant values above
#' confidence limits of random expectations. Default = "#000000".
#' @param col_random_values color code or name for randomized values.
#' Default = "D2D2D2".
#' @param pch_all point symbol to be used for all values. Defaults = 1 or 16.
#' @param pch_significant_low point symbol to be used for significant values
#' below confidence limits of random expectations. Default = 16.
#' @param pch_significant_high point symbol to be used for significant values
#' above confidence limits of random expectations. Default = 16.
#' @param pch_random_values point symbol to be used for randomized values.
#' Default = 1.
#' @param main main title for the plot. Default = NULL.
#' @param xlab label for the x-axis. Default = NULL.
#' @param ylab label for the y-axis. Default = NULL.
#' @param xlim x limits of the plot (x1, x2). For \code{plot_PAM_CS}, the
#' default, NULL, uses the range of normalized richness.
#' @param ylim y limits of the plot. For \code{plot_PAM_CS}, the
#' default, NULL, uses the range of the normalized values of the dispersion
#' field. The second limit is increased by adding the result of multiplying it
#' by \code{ylim_expansion}, if \code{add_legend} = TRUE.
#' @param ylim_expansion value used or expanding the \code{ylim}. Default =
#' 0.25.
#' @param las the style of axis labels; default = 1.
#' See \code{\link[graphics]{par}}.
#' @param add_legend (logical) whether to add a legend describing information
#' relevant for interpreting the diagram. Default = TRUE.
#'
#' @return
#' For \code{plot_PAM_CS}:
#'
#' A range-diversity plot with values of normalized richness in the x-axis,
#' and normalized values of the dispersion field index divided by number of
#' species in the y-axis.
#'
#' For \code{plot_PAM_CS_geo}:
#'
#' A geographic view of the PAM representing the areas or points identified
#' as non statistically significant, significant above random expectations,
#' and significant below random expectations.
#'
#' @usage
#' plot_PAM_CS(PAM_CS, add_significant = FALSE,
#'             add_random_values = FALSE, col_all = "#CACACA",
#'             col_significant_low = "#6D6D6D",
#'             col_significant_high = "#000000",
#'             col_random_values = "#D2D2D2", pch_all = 1,
#'             pch_significant_low = 16, pch_significant_high = 16,
#'             pch_random_values = 1, main = NULL,
#'             xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
#'             ylim_expansion = 0.25, las = 1, add_legend = TRUE)
#'
#' @export
#' @importFrom graphics legend polygon
#'
#' @rdname plot_PAM_CS
#'
#' @examples
#' # Data
#' b_pam <- read_PAM(system.file("extdata/b_pam.rds",
#'                               package = "biosurvey"))
#'
#' # Preparing data for CS diagram
#' pcs <- prepare_PAM_CS(PAM = b_pam)
#'
#' # Plot
#' plot_PAM_CS(pcs)

plot_PAM_CS <- function(PAM_CS, add_significant = FALSE,
                        add_random_values = FALSE, col_all = "#CACACA",
                        col_significant_low = "#6D6D6D",
                        col_significant_high = "#000000",
                        col_random_values = "#D2D2D2", pch_all = 1,
                        pch_significant_low = 16, pch_significant_high = 16,
                        pch_random_values = 1, main = NULL,
                        xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
                        ylim_expansion = 0.25, las = 1, add_legend = TRUE) {

  # Initial tests
  if (!class(PAM_CS)[1] %in% c("base_PAM", "PAM_CS")) {
    stop("Class of 'PAM_CS' is not supported, see function's help.")
  }

  # Preparing data
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

  # Prepare vertex of plot limits
  vx <- PAM$Theoretical_boundaries$x
  vy <- PAM$Theoretical_boundaries$y
  sper <- round(PAM$Spearman_cor, 3)

  # Plot elements
  if (is.null(main)) {
    main <- "Range-diversity plot"
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

  # Plot
  plot(alfas, fists, col = col_all, pch = pch_all, xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab, main = main, las = las)
  polygon(vx, vy, border = "#474747")
  if (add_legend == TRUE) {
    legend("topleft", bty = "n", inset = -0.02,
           legend = c(paste("N species (S) =", s), paste("N sites-cells =", n),
                      as.expression(bquote("Whittaker's" ~ beta ~ "=" ~
                                             .(round(betty, 3)))),
                      as.expression(bquote("Spearman's" ~ r[s] ~ "=" ~
                                             .(sper)))))
  }

  # Randomized values
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

  # Significant values
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




#' @rdname plot_PAM_CS
#' @export
#' @param xy_coordinates Only required if \code{PAM_CS} is of class PAM_CS.
#' A matrix or data.frame containing the columns longitude and latitude (in that
#' order) corresponding to the points in \code{PAM_CS$S_significance_id}.
#' Default = NULL
#' @param border color for cell borders of the PAM grid. The default, NULL, does
#' not plot any border.
#' @param mar (numeric) vector of length 4 to set the margins of the plot in
#' geography. The default, NULL, is (3.1, 3.1, 2.1, 2.1).
#'
#' @usage
#' plot_PAM_CS_geo(PAM_CS, xy_coordinates = NULL, col_all = "#CACACA",
#'                 col_significant_low = "#6D6D6D",
#'                 col_significant_high = "#000000", border = NULL,
#'                 pch_all = 16, pch_significant_low = 16,
#'                 pch_significant_high = 16, xlim = NULL,
#'                 ylim = NULL, mar = NULL)

plot_PAM_CS_geo <- function(PAM_CS, xy_coordinates = NULL, col_all = "#CACACA",
                            col_significant_low = "#6D6D6D",
                            col_significant_high = "#000000", border = NULL,
                            pch_all = 16, pch_significant_low = 16,
                            pch_significant_high = 16, xlim = NULL,
                            ylim = NULL, mar = NULL) {

  # Initial tests
  if (!class(PAM_CS)[1] %in% c("base_PAM", "PAM_CS")) {
    stop("Class of 'PAM_CS' is not supported, see function's help.")
  }

  # Preparing data
  if (class(PAM_CS)[1] == "base_PAM") {
    bp <- PAM_CS$PAM
    PAM <- PAM_CS$PAM_indices$CS_diagram

    # additional tests
    if (all(is.na(PAM$S_significance_id))) {
      stop("Values that indicate significance are missing in 'PAM_CS'")
    }

    ## Box to plot
    boxp <- matrix(terra::ext(bp), nrow = 2)
    boxpam <- terra::vect(boxp, crs = terra::crs(bp))

    ## limits
    if (is.null(xlim)) {
      xlim <- boxp[, 1]
    }
    if (is.null(ylim)) {
      ylim <- boxp[, 2]
    }
  } else {
    PAM <- PAM_CS

    # additional tests
    if (all(is.na(PAM$S_significance_id))) {
      stop("Values that indicate significance are missing in 'PAM_CS'")
    }

    if (is.null(xy_coordinates)) {
      stop("If argument 'PAM_CS' is of class 'PAM_CS', 'xy_coordinates' must be defined")
    } else {
      if (!class(xy_coordinates)[1] %in% c("matrix", "data.frame")) {
        stop("'xy_coordinates' must be of class 'matrix' or 'data.frame'")
      }

      if (nrow(xy_coordinates) != length(PAM$S_significance_id)) {
        stop("Number of 'xy_coordinates' and points in 'PAM_CS' do not match")
      }
    }

    # limits
    if (is.null(xlim)) {
      xlim <- range(xy_coordinates[, 1])
      dif <- diff(xlim) * 0.05
      xlim <- xlim + c(-dif, dif)
    }
    if (is.null(ylim)) {
      ylim <- range(xy_coordinates[, 2])
      dif <- diff(ylim) * 0.05
      ylim <- ylim + c(-dif, dif)
    }
  }

  # Colors
  sigfact <- as.factor(PAM$S_significance_id)
  cols <- c(col_all, col_significant_low, col_significant_high)

  # Plot
  if (class(PAM_CS)[1] == "base_PAM") {
    ## bp grid
    bp <- bp[bp$ID %in% names(PAM$S_significance_id), ]

    ## border
    border <- ifelse(is.null(border), NA, border)

    ## The plot
    terra::plot(boxpam, col = NA, xlim = xlim, ylim = ylim,
                axes = FALSE, legend = FALSE, mar = mar)
    maps::map(fill = TRUE, col = "gray98", lforce = "n",
              border = "gray90", add = TRUE)
    terra::plot(bp, col = cols[sigfact], border = border, add = TRUE,
                legend = FALSE)
    box()
  } else {
    pchs <- c(pch_all, pch_significant_low, pch_significant_high)

    maps::map(fill = TRUE, col = "gray98", lforce = "n",
              border = "gray90", xlim = xlim, ylim = ylim)
    box()
    points(xy_coordinates, col = cols[sigfact], pch = pchs[sigfact])
  }
}
