#' Representation of environmental blocks in geography and environment
#'
#' @description creates a two-panel plot representing environmental blocks
#' (all or selected) in both spaces, environmental and geographic.
#'
#' @param master a master_matrix object derived from function
#' \code{\link{prepare_master_matrix}} or a master_selection object derived
#' from functions \code{\link{uniformG_selection}},
#' \code{\link{uniformE_selection}} or \code{\link{EG_selection}}.
#' @param variable_1 (character or numeric) name or position of the first
#' variable (X-axis) used to create blocks.
#' @param variable_2 (character or numeric) name or position of the second
#' variable (Y-axis) used to create blocks (must be different from the
#' first one).
#' @param region_border (logical) whether to add region border to the plot.
#' Default = TRUE.
#' @param mask_border (logical) whether to add mask border to the plot. Ignored
#' if mask is not present in \code{master_selection}. Default = FALSE.
#' @param which (character) blocks to be plotted. Options are "all" or
#' "selected". Default = "all".
#' @param block_ID (logical) whether to add a text ID to blocks plotted in
#' environmental space. Default = FALSE.
#' @param col_all colors for points in all blocks. The default, NULL, uses a
#' color blind friendly palette to differentiate among distinct blocks when
#' \code{which} = "all", or uses a light gray color when
#' \code{which} = "selected". See details for explanations of how to define
#' them.
#' @param col_selected color for points in selected blocks. Ignored if
#' \code{which} = "all". The default, NULL, uses a blue color to represent
#' selected blocks on top of all.
#' @param col_ID color for text ID to be added if \code{block_ID} = TRUE. The
#' default, NULL, uses the "back".
#' @param cex_all (numeric) value defining magnification of points in all blocks
#' relative to the default. Default = 0.7.
#' @param cex_selected (numeric) value defining magnification of points in
#' selected blocks relative to the default. Default = 1.
#' @param cex_ID (numeric) value defining magnification of text ID to be added
#' if \code{block_ID} = TRUE. Default = 1.
#' @param pch_all (numeric) integer specifying a symbol when plotting points of
#' all blocks. Default = 16.
#' @param pch_selected (numeric) integer specifying a symbol when plotting
#' points of selected blocks. Default = 16.
#' @param add_main (logical) whether or not to add fixed titles to the plot.
#' Default = TRUE. Titles added are "Environmental space" and "Geographic
#' space".
#'
#' @details
#' Defining colors in \code{col_all} depends on what is chosen in \code{which}.
#' If "all" is chosen, it is convenient to define \code{col_all} as a color
#' ramp palette (randomly arranged) or a set of colors depending on the number
#' of blocks in the object defined in \code{master}. If "selected" is chosen
#' in \code{which} it is recommended to use a single color, preferably a light
#' one, so the selected blocks can be easily identified. See examples.
#'
#' @return
#' A two-panel plot showing all the blocks of the region of interest and the
#' blocks that were selected. They are show in both spaces, geographic and
#' environmental.
#'
#' @usage
#' plot_blocks_EG(master, variable_1, variable_2, region_border = TRUE,
#'                mask_border = FALSE, which = "all", block_ID = FALSE,
#'                col_all = NULL, col_selected = NULL, col_ID = NULL,
#'                cex_all = 0.7, cex_selected = 1, cex_ID = 1,
#'                pch_all = 16, pch_selected = 16, add_main = TRUE)
#'
#' @export
#' @importFrom maps map
#' @importFrom sp plot
#' @importFrom graphics layout par plot.new text title
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#'
#' # Creating blocks
#' m_blocks <- make_blocks(m_matrix, variable_1 = "PC1",
#'                         variable_2 = "PC2", n_cols = 10, n_rows = 10,
#'                         block_type = "equal_area")
#'
#' plot_blocks_EG(master = m_blocks, variable_1 = "PC1", variable_2 = "PC2",
#'                block_ID = TRUE)
#'
#' # Defining your own colors
#' n_blocks <- length(m_blocks$data_matrix$Block)
#' your_palette <- sample(heat.colors(n_blocks), n_blocks)
#' block_factor <- as.factor(m_blocks$data_matrix$Block)
#' your_colors <- your_palette[block_factor]
#'
#' plot_blocks_EG(master = m_blocks, variable_1 = "PC1", variable_2 = "PC2",
#'                block_ID = TRUE, col_all = your_colors)


plot_blocks_EG <- function(master, variable_1, variable_2, region_border = TRUE,
                           mask_border = FALSE,  which = "all", block_ID = FALSE,
                           col_all = NULL, col_selected = NULL, col_ID = NULL,
                           cex_all = 0.7, cex_selected = 1, cex_ID = 1,
                           pch_all = 16, pch_selected = 16, add_main = TRUE) {
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
  if (!which %in% c("all", "selected")) {
    stop("Argument 'which' is not valid, options are 'all' or 'selected'.")
  } else {
    if (is.null(master$data_matrix$Block)) {
      stop("Blocks are not defined in data_matrix of 'master', see function 'make_blocks'.")
    }
    if (which == "selected") {
      if (is.null(master$data_matrix$Selected_blocks)) {
        stop("Object in 'master' does not contain selected blocks, see function 'block_sample'.")
      }
    }
  }

  # preparing data
  ublocks <- unique(master$data_matrix$Block)
  nblocks <- length(ublocks)

  gvars <- c("Longitude", "Latitude")
  evars <- c(variable_1, variable_2)

  # colors
  col_pal <- colorRampPalette(rev(c("#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb",
                                    "#41b6c4", "#1d91c0", "#225ea8", "#253494",
                                    "#081d58")))

  if (is.null(col_all) & is.null(col_selected) & is.null(col_ID)) {
    if (which == "all") {
      col_all <- sample(col_pal(nblocks),
                        nblocks)[as.factor(master$data_matrix$Block)]
    } else {
      col_all <- "#E1E1E1"
      col_selected <- "#3B22CB"
    }
    col_ID <- "#000000"
  } else {
    if (is.null(col_all)) {
      if (which == "all") {
        col_all <- sample(col_pal(nblocks), nblocks)
      } else {
        col_all <- "#E1E1E1"
      }
    }
    if (is.null(col_selected)) {
      col_selected <- "#3B22CB"
    }
    if (is.null(col_ID)) {
      col_ID <- "#000000"
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
  plot(master$data_matrix[, evars], col = col_all, pch = pch_all, cex = cex_all,
       bty = "l", xlab = "", ylab = "")
  title(xlab = variable_1, line = 2.3, cex.lab = 1.1)
  title(ylab = variable_2, line = 2.3, cex.lab = 1.1)

  ## selected blocks
  if (which == "selected") {
    sel <- which(master$data_matrix$Selected_blocks == 1)
    selected_data <- master$data_matrix[sel, ]
    points(selected_data[sel, evars], pch = pch_selected, cex = cex_selected,
           col = col_selected)
  }

  ## add block ID
  if (block_ID == TRUE) {
    cents <- lapply(ublocks, function(x) {
      cen <- apply(master$data_matrix[master$data_matrix$Block == x, evars],
                   2, mean)
      text(cen[1], cen[2], labels = x, cex = cex_ID, col = col_ID)
    })
  }

  ## geographic space
  ### box to plot
  boxpam <- t(master$region@bbox)
  boxpam <- sp::SpatialPointsDataFrame(boxpam, data.frame(boxpam),
                                       proj4string = master$region@proj4string)

  ### plot
  par(mar = rep(0.5, 4))
  sp::plot(boxpam, col = NA)
  maps::map(fill = TRUE, col = "gray97", lforce = "n",
            border = "gray80", add = TRUE)
  box(which = "plot")
  points(master$data_matrix[, gvars], pch = pch_all, cex = cex_all,
         col = col_all)

  if (region_border == TRUE) {
    sp::plot(master$region, border = "gray50", add = TRUE)
  }
  if (mask_border == TRUE & !is.null(master$mask)) {
    sp::plot(master$mask, border = "gray50", add = TRUE)
  }

  ## selected blocks
  if (which == "selected") {
    points(selected_data[sel, gvars], pch = pch_selected, cex = cex_selected,
           col = col_selected)
  }
}

