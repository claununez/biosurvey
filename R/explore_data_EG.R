#' Plots to explore environmental factors in environmental and geographic space
#'
#' @description Creates a four-panel plot with information of two environmental
#' predictors (at a time) in the region of interest. The two top panels contain
#' the information in geographic space (one predictor per panel). The two panels
#' at the bottom contain information in a 2D environmental space for the two
#' variables.
#'
#' @param master a master_matrix object derived from function
#' \code{\link{master_matrix}} or a master_selection object derived from functions
#' \code{\link{uniformG_selection}}, \code{\link{uniformE_selection}}
#' or \code{EG_selection}.
#' @param variable_1 (character or numeric) name or position of the first
#' variable (X axis) to be explored.
#' @param variable_2 (character or numeric) name or position of the second
#' variable (Y axis) to be explored (must be different from the first one).
#' @param col_variable1 a color palette for \code{variable_1} defined
#' using functions like \code{\link[grDevices]{heat.colors}}, or one generated
#' using functions like \code{\link[grDevices]{colorRampPalette}}. The default,
#' NULL, uses \code{viridis::cividis(255)}.
#' @param col_variable2 a color palette for \code{variable_2} defined
#' using functions like \code{\link[grDevices]{heat.colors}}, or one generated
#' using functions like \code{\link[grDevices]{colorRampPalette}}. The default,
#' NULL, uses \code{viridis::cividis(255)}.
#' @param col_points color for points in environmental space. The default, NULL,
#' uses the 25th color of the default palette for \code{col_variable1} with an
#' alpha of 0.6.
#' @param col_density a color palette to represent representation density of points
#' in environmental space. This palette can be defined using functions like
#' \code{\link[grDevices]{heat.colors}}, or one generated using functions like
#' \code{\link[grDevices]{colorRampPalette}}. The default, NULL, uses
#' \code{viridis::viridis(255)}, and changes the first color in the palette by NA.
#'
#' @return
#' A multipanel plot showing two of the environmental preditors in the region of
#' interest in both spaces, geographic and environmental.
#'
#' @usage
#' explore_data_EG(master, variable_1, variable_2, col_variable1 = NULL,
#'                 col_variable2 = NULL, col_points = NULL, col_density = NULL)
#'
#' @export
#' @importFrom ks kde
#' @importFrom viridis viridis cividis
#' @importFrom scales alpha
#' @importFrom sp plot
#' @importFrom raster image
#' @importFrom graphics image layout par plot.new rasterImage text title
#' @importFrom grDevices as.raster
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#'
#' colnames(m_matrix$master_matrix)
#'
#' # Plot
#' explore_data_EG(m_matrix, variable_1 = "Mean_temperature",
#'                 variable_2 = "Annual_precipitation")


explore_data_EG <- function(master, variable_1, variable_2,
                            col_variable1 = NULL, col_variable2 = NULL,
                            col_points = NULL, col_density = NULL) {

  # initial tests
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

  # par settings
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  # kernel
  mx2kd <- ks::kde(master$master_matrix[, c(variable_1, variable_2)])

  # colors
  if (is.null(col_variable1) & is.null(col_variable2) & is.null(col_points) &
      is.null(col_density)) {
    col_density <- viridis::viridis(255)
    col_points <- scales::alpha(col_density[25], 0.6)
    col_variable1 <- viridis::cividis(255)
    col_variable2 <- col_variable1
    col_density[1] <- NA
  } else {
    if (is.null(col_variable1)) {
      col_variable1 <- viridis::cividis(255)
    }
    if (is.null(col_variable2)) {
      col_variable2 <- viridis::cividis(255)
    }
    if (is.null(col_variable2)) {
      col_variable2 <- viridis::cividis(255)
    }
    if (is.null(col_density)) {
      col_density <- viridis::viridis(255)
      col_density[1] <- NA
    }
    if (is.null(col_points)) {
      col_points <- scales::alpha(viridis::viridis(255)[25], 0.6)
    }
  }

  # limits
  xlim <- range(master$master_matrix[, variable_1])
  ylim  <- range(master$master_matrix[, variable_2])

  # plot
  layout(matrix(1:20, 4, byrow = T), widths = c(1, 10, 2, 10, 2),
                   heights = c(1, 10, 1.2, 10))
  par(cex = 0.7, mar = rep(0, 4))

  ## titles
  plot.new()
  plot.new()
  text(0.5, 0.5, variable_1, cex = 1.2)
  plot.new()
  plot.new()
  text(0.5, 0.5, variable_2, cex = 1.2)
  plot.new()

  ## geographic space
  plot.new()
  text(0.5, 0.5, "Geographic space", cex = 1.2, srt = 90)

  sp::plot(master$polygon, border = NA)
  var <- master$raster_base
  var[!is.na(var[])] <- master$master_matrix[, variable_1]
  value_range <- c(var@data@min, var@data@max)
  raster::image(var, col = col_variable1, add = TRUE)
  sp::plot(master$polygon, add = TRUE)

  plot.new()
  bar_legend(value_range, col = col_variable1, title = variable_1)


  sp::plot(master$polygon, border = NA)
  var[!is.na(var[])] <- master$master_matrix[, variable_2]
  value_range <- c(var@data@min, var@data@max)
  raster::image(var, col = col_variable2, add = TRUE)
  sp::plot(master$polygon, add = TRUE)

  plot.new()
  bar_legend(value_range, col = col_variable2, title = variable_2)

  ## titles
  plot.new()
  plot.new()
  text(0.5, 0.5, "Point cloud", cex = 1.2)
  plot.new()
  plot.new()
  text(0.5, 0.5, "Point density", cex = 1.2)
  plot.new()

  ## environmental space
  plot.new()
  text(0.5, 0.6, "Environmental space", cex = 1.2, srt = 90)

  par(mar = c(3.5, 3.5, 0.5, 0.5))
  plot(master$master_matrix[, c(variable_1, variable_2)], col = col_points,
       bty = "l", xlab = "", ylab = "")
  title(xlab = variable_1, line = 2.4, cex.lab = 1.1)
  title(ylab = variable_2, line = 2.4, cex.lab = 1.1)

  par(mar = rep(0, 4))
  plot.new()

  par(mar = c(3.5, 3.5, 0.5, 0.5))
  plot(xlim, ylim, type = "n", bty = "l", xlab = "", ylab = "")
  image(mx2kd$eval.points[[1]], mx2kd$eval.points[[2]],
                  z = mx2kd$estimate, col = col_density, add = TRUE)
  title(xlab = variable_1, line = 2.4, cex.lab = 1.1)

  par(mar = rep(0, 4))
  plot.new()
  bar_legend(c("Low", "High"), col = col_density, title = "Density", round = 1)
}


