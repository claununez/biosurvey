#' Plots to explore environmental factors in environmental and geographic space
#'
#' @description Creates a four-panel plot with information of two environmental
#' predictors (at a time) in the region of interest (or region reduced with mask,
#' if used). The two top panels contain the information in geographic space (one
#' predictor per panel). The two panels at the bottom contain information in a
#' 2D environmental space for the two variables.
#'
#' @param master master_matrix object derived from function
#' \code{\link{prepare_master_matrix}} or master_selection object derived from
#' functions \code{\link{uniformG_selection}}, \code{\link{uniformE_selection}}
#' or \code{\link{EG_selection}}.
#' @param variable_1 (character or numeric) name or position of the first
#' variable (x-axis) to be explored.
#' @param variable_2 (character or numeric) name or position of the second
#' variable (y-axis) to be explored (must be different from the first one).
#' @param region_border (logical) whether to add region border to the plot.
#' Default = TRUE.
#' @param mask_border (logical) whether to add mask border to the plot. Ignored
#' if mask is not present in \code{master_selection}. Default = FALSE.
#' @param col_variable1 a color palette for \code{variable_1} defined
#' using functions like \code{\link{heat.colors}}, or one generated
#' using functions like \code{\link{colorRampPalette}}. The default,
#' NULL, uses a colorblind friendly palette similar to viridis.
#' @param col_variable2 a color palette for \code{variable_2} defined
#' using functions like \code{\link{heat.colors}}, or one generated
#' using functions like \code{\link{colorRampPalette}}. The default,
#' NULL, uses a colorblind friendly palette similar to viridis.
#' @param col_points color for points in environmental space. The default, NULL,
#' uses the 25th color of the default palette for \code{col_variable1} with an
#' alpha of 0.6.
#' @param col_density color palette to represent representation density of
#' points in environmental space. This palette can be defined using functions
#' like \code{\link{heat.colors}}, or one generated using functions like
#' \code{\link{colorRampPalette}}. The default, NULL, uses a colorblind
#' friendly palette similar to magma, and changes the first color in the palette
#' to NA.
#'
#' @return
#' A multi-panel plot showing two of the environmental predictors in the region
#' of interest in both spaces, geographic and environmental.
#'
#' @usage
#' explore_data_EG(master, variable_1, variable_2, region_border = TRUE,
#'                 mask_border = FALSE, col_variable1 = NULL,
#'                 col_variable2 = NULL, col_points = NULL, col_density = NULL)
#'
#' @export
#' @importFrom ks kde
#' @importFrom sp plot
#' @importFrom raster image
#' @importFrom graphics image layout par plot.new rasterImage text title box points
#' @importFrom grDevices as.raster col2rgb rgb
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#'
#' colnames(m_matrix$data_matrix)
#'
#' # Plot
#' \donttest{
#' explore_data_EG(m_matrix, variable_1 = "Mean_temperature",
#'                 variable_2 = "Annual_precipitation")
#' }


explore_data_EG <- function(master, variable_1, variable_2,
                            region_border = TRUE, mask_border = FALSE,
                            col_variable1 = NULL, col_variable2 = NULL,
                            col_points = NULL, col_density = NULL) {

  # Initial tests
  if (missing(master)) {
    stop("Argument 'master' is required to produce the plot.")
  }
  if (missing(variable_1)) {
    stop("Argument 'variable_1' is required to produce the plot.")
  }
  if (missing(variable_2)) {
    stop("Argument 'variable_2' is required to produce the plot.")
  }
  if (!class(master)[1] %in% c("master_matrix", "master_selection")) {
    stop("Object defined in 'master' is not valid, see function's help.")
  }

  # Par settings
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  # Kernel
  mx2kd <- ks::kde(master$data_matrix[, c(variable_1, variable_2)])

  # Colors
  ## Palettes
  col_pal1 <- purplow
  col_pal2 <- darkros

  if (is.null(col_variable1) & is.null(col_variable2) & is.null(col_points) &
      is.null(col_density)) {
    ## If all null
    col_density <- col_pal2(255)
    col_points <- make_alpha(col_density[25], 0.6)
    uv1 <- length(unique(master$data_matrix[, variable_1]))
    col_variable1 <- col_pal1(uv1)
    uv2 <- length(unique(master$data_matrix[, variable_2]))
    col_variable2 <- col_pal1(uv2)
    col_density[1] <- NA
  } else {
    ## If some of them are null
    if (is.null(col_variable1)) {
      uv1 <- length(unique(master$data_matrix[, variable_1]))
      col_variable1 <- col_pal1(uv1)
    }
    if (is.null(col_variable2)) {
      uv2 <- length(unique(master$data_matrix[, variable_2]))
      col_variable2 <- col_pal1(uv2)
    }
    if (is.null(col_density)) {
      col_density <- col_pal2(255)
      col_density[1] <- NA
    }
    if (is.null(col_points)) {
      col_points <- make_alpha(col_pal2(255)[25], 0.6)
    }
  }

  # Limits
  xlim <- range(master$data_matrix[, variable_1])
  ylim  <- range(master$data_matrix[, variable_2])

  # Box to plot
  boxpam <- t(master$region@bbox)
  boxpam <- sp::SpatialPointsDataFrame(boxpam, data.frame(boxpam),
                                       proj4string = master$region@proj4string)

  # to fill cells
  to_fill <- !is.na(master$raster_base[])

  # Plot
  layout(matrix(1:15, 3, byrow = T), widths = c(1, 10, 2, 10, 2),
                   heights = c(10, 1.2, 10))
  par(cex = 0.7, mar = rep(0, 4))

  ## Geographic space
  plot.new()
  text(0.5, 0.5, "Geographic space", cex = 1.2, srt = 90)

  ### Variable 1
  #vf <- as.factor(master$data_matrix[, variable_1])
  #sp::plot(master$raster_base, col = col_variable1[vf], border = NA)
  master$raster_base[to_fill] <- master$data_matrix[, variable_1]
  sp::plot(boxpam, col = NA)
  raster::image(master$raster_base, col = col_variable1, add = TRUE)

  if(region_border == TRUE) {
    sp::plot(master$region, border = "gray50", add = TRUE)
  }
  if (mask_border == TRUE & !is.null(master$mask)) {
    sp::plot(master$mask, border = "gray50", add = TRUE)
  }

  plot.new()
  bar_legend(xlim, col = col_variable1, title = variable_1)

  ### Variable 2
  #vf <- as.factor(master$data_matrix[, variable_2])
  #sp::plot(master$raster_base, col = col_variable2[vf], border = NA)
  master$raster_base[to_fill] <- master$data_matrix[, variable_2]
  sp::plot(boxpam, col = NA)
  raster::image(master$raster_base, col = col_variable2, add = TRUE)

  if(region_border == TRUE) {
    sp::plot(master$region, border = "gray50", add = TRUE)
  }
  if (mask_border == TRUE & !is.null(master$mask)) {
    sp::plot(master$mask, border = "gray50", add = TRUE)
  }

  plot.new()
  bar_legend(ylim, col = col_variable2, title = variable_2)

  ## Titles
  plot.new()
  plot.new()
  text(0.55, 0.5, "Point cloud", cex = 1.2)
  plot.new()
  plot.new()
  text(0.55, 0.5, "Point density", cex = 1.2)
  plot.new()

  ## Environmental space
  plot.new()
  text(0.5, 0.6, "Environmental space", cex = 1.2, srt = 90)

  par(mar = c(3.5, 3.5, 0.5, 0.5))
  plot(master$data_matrix[, c(variable_1, variable_2)], col = col_points,
       bty = "l", xlab = "", ylab = "")
  title(xlab = variable_1, line = 2.1, cex.lab = 1.1)
  title(ylab = variable_2, line = 2.1, cex.lab = 1.1)

  par(mar = rep(0, 4))
  plot.new()

  par(mar = c(3.5, 3.5, 0.5, 0.5))
  plot(xlim, ylim, type = "n", bty = "l", xlab = "", ylab = "")
  image(mx2kd$eval.points[[1]], mx2kd$eval.points[[2]],
                  z = mx2kd$estimate, col = col_density, add = TRUE)
  title(xlab = variable_1, line = 2.1, cex.lab = 1.1)

  par(mar = rep(0, 4))
  plot.new()
  bar_legend(c("Low", "High"), col = col_density, title = "Density", round = 1)
}


