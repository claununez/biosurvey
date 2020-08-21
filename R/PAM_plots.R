#' Plot of PAM indices in geography
#'
#' @param PAM an object of class base_PAM.
#' @param index (character) code for the index to be plotted. Options are: "RI"
#' (Richness), "RIN" (Richness normalized), "DF" (Dispersion field), or "MCC"
#' (Mean composition covariance). Default = "RI".
#' @param col_pal color palette function to be used in defining colors for the
#' \code{index} to be plotted. The default, NULL, uses a color blind friendly
#' palette similar to viridis.
#' @param border color for cell borders of the PAM grid. The default, NULL, does
#' not plot any border.
#' @param cex (numeric) value by which plotting elements should be magnified
#' relative to the default. Default = 0.9
#'
#' @return
#' A plot of \code{index} represented in geography.
#'
#' @usage
#' plot_PAM_geo(PAM, index = "RI", col_pal = NULL, border = NULL, cex = 0.9)
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

plot_PAM_geo <- function(PAM, index = "RI", col_pal = NULL, border = NULL,
                         cex = 0.9) {
  if (missing(PAM)) {
    stop("Argument 'PAM' is missing.")
  }
  all_in <- c("RI", "RIN", "DF", "MCC")
  if (!index %in% all_in) {
    stop("Argument 'index' is not valid, options are: 'RI', 'RIN', 'DF', or 'MCC'.")
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
  plot.new()
  bar_legend(rfactor, col = col, title = gsub("_", " ", g_indices[index]),
             round = 3, label_x = 0.5, labels_y = c(0.18, 0.87))
}



