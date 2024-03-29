#' Representation of environmental blocks in geography and environment
#'
#' @description Creates a plot representing environmental blocks
#' (all or selected) in both spaces, environmental and/or geographic.
#'
#' @param master master_matrix object derived from function
#' \code{\link{prepare_master_matrix}} or master_selection object derived
#' from functions \code{\link{uniformG_selection}},
#' \code{\link{uniformE_selection}} or \code{\link{EG_selection}}. Blocks must
#' be defined, see \code{\link{make_blocks}}.
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
#' @param mar (numeric) vector of length 4 to set the margins of the plot in G.
#' The default, NULL, is (3.1, 3.1, 2.1, 2.1) for `plot_blocks_G` and
#' (3.5, 0.5, 0.5, 0.5) for `plot_blocks_EG.`
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
#' A plot showing all the blocks of the region of interest and, if asked, the
#' blocks that were selected. They are show in both spaces, geographic and/or
#' environmental.
#'
#' @usage
#' plot_blocks_EG(master, region_border = TRUE, mask_border = FALSE,
#'                which = "all", block_ID = FALSE, col_all = NULL,
#'                col_selected = NULL, col_ID = NULL, cex_all = 0.7,
#'                cex_selected = 1, cex_ID = 1, pch_all = 16,
#'                pch_selected = 16, add_main = TRUE, mar = NULL)
#'
#' @export
#' @importFrom maps map
#' @importFrom terra plot ext vect crs
#' @importFrom graphics layout par plot.new text title
#'
#' @examples
#' # Data
#' m_matrix <- read_master(system.file("extdata/m_matrix.rds",
#'                                     package = "biosurvey"))
#'
#' # Creating blocks
#' m_blocks <- make_blocks(m_matrix, variable_1 = "PC1",
#'                         variable_2 = "PC2", n_cols = 10, n_rows = 10,
#'                         block_type = "equal_area")
#'
#' plot_blocks_EG(master = m_blocks, block_ID = TRUE)
#' plot_blocks_E(master = m_blocks)
#' plot_blocks_G(master = m_blocks)
#'
#' # Defining your own colors
#' n_blocks <- length(m_blocks$data_matrix$Block)
#' your_palette <- sample(heat.colors(n_blocks), n_blocks)
#' block_factor <- as.factor(m_blocks$data_matrix$Block)
#' your_colors <- your_palette[block_factor]
#'
#' plot_blocks_EG(master = m_blocks, block_ID = TRUE, col_all = your_colors)


plot_blocks_EG <- function(master, region_border = TRUE, mask_border = FALSE,
                           which = "all", block_ID = FALSE,
                           col_all = NULL, col_selected = NULL, col_ID = NULL,
                           cex_all = 0.7, cex_selected = 1, cex_ID = 1,
                           pch_all = 16, pch_selected = 16, add_main = TRUE,
                           mar = NULL) {
  # Initial tests
  if (missing(master)) {
    stop("Argument 'master' is required to produce the plot.")
  }
  if (!class(master)[1] %in% c("master_matrix", "master_selection")) {
    stop("Object defined in 'master' is not valid, see function's help.")
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
  plot_blocks_E(master, which, block_ID, col_all,
                col_selected, col_ID, cex_all, cex_selected, cex_ID,
                pch_all, pch_selected)

  ## Geographic space
  if (is.null(mar)) {
    mar <- c(3.5, rep(0.5, 3))
  }
  plot_blocks_G(master, region_border, mask_border, which, block_ID, col_all,
                col_selected, col_ID, cex_all, cex_selected, cex_ID, pch_all,
                pch_selected, mar)
}


#' @rdname plot_blocks_EG
#' @export
#' @param main (character) the main title for the plot.
#' @param xlab (character) label for the x axis. The default, NULL, uses
#' variable_1.
#' @param ylab (character) label for the y axis. The default, NULL, uses
#' variable_2.
#' @usage
#' plot_blocks_E(master, which = "all", block_ID = FALSE, col_all = NULL,
#'               col_selected = NULL, col_ID = NULL, cex_all = 0.7,
#'               cex_selected = 1, cex_ID = 1, pch_all = 16,
#'               pch_selected = 16, main = "", xlab = NULL, ylab = NULL)

plot_blocks_E <- function(master, which = "all", block_ID = FALSE,
                          col_all = NULL, col_selected = NULL,
                          col_ID = NULL, cex_all = 0.7, cex_selected = 1,
                          cex_ID = 1, pch_all = 16, pch_selected = 16,
                          main = "", xlab = NULL, ylab = NULL) {
  # Initial tests
  if (missing(master)) {
    stop("Argument 'master' is required to produce the plot.")
  }
  if (!class(master)[1] %in% c("master_matrix", "master_selection")) {
    stop("Object defined in 'master' is not valid, see function's help.")
  }
  if (is.null(master$data_matrix$Block)) {
    stop("Blocks are not defined in data_matrix, see function 'make_blocks'.")
  } else {
    sel_args <- attributes(master$data_matrix)

    variable_1 <- sel_args$arguments$variable_1
    variable_2 <- sel_args$arguments$variable_2

    coln <- colnames(master$data_matrix)
    if (!variable_1 %in% coln) {
      stop(variable_1, " is not one o the columns in 'master$data_matrix'.")
    }
    if (!variable_2 %in% coln) {
      stop(variable_2, " is not one o the columns in 'master$data_matrix'.")
    }
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

  # Preparing data
  ublocks <- unique(master$data_matrix$Block)
  nblocks <- length(ublocks)

  evars <- c(variable_1, variable_2)

  # Colors
  col_pal <- purplow

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

  # labels
  if (is.null(xlab)) {
    xlab <- variable_1
  }
  if (is.null(ylab)) {
    ylab <- variable_2
  }

  ## Plot
  plot(master$data_matrix[, evars], col = col_all, pch = pch_all, cex = cex_all,
       bty = "l", xlab = "", ylab = "", main = main)
  title(xlab = xlab, line = 2.3, cex.lab = 1.1)
  title(ylab = ylab, line = 2.3, cex.lab = 1.1)

  ## Selected blocks
  if (which == "selected") {
    sel <- which(master$data_matrix$Selected_blocks == 1)
    selected_data <- master$data_matrix[sel, ]
    points(selected_data[sel, evars], pch = pch_selected, cex = cex_selected,
           col = col_selected)
  }

  ## Add block ID
  if (block_ID == TRUE) {
    cents <- lapply(ublocks, function(x) {
      cen <- apply(master$data_matrix[master$data_matrix$Block == x, evars],
                   2, mean)
      text(cen[1], cen[2], labels = x, cex = cex_ID, col = col_ID)
    })
  }
}



#' @rdname plot_blocks_EG
#' @export
#' @usage
#' plot_blocks_G(master, region_border = TRUE, mask_border = FALSE,
#'               which = "all", block_ID = FALSE, col_all = NULL,
#'               col_selected = NULL, col_ID = NULL, cex_all = 0.7,
#'               cex_selected = 1, cex_ID = 1, pch_all = 16,
#'               pch_selected = 16, mar = NULL)

plot_blocks_G <- function(master, region_border = TRUE, mask_border = FALSE,
                          which = "all", block_ID = FALSE, col_all = NULL,
                          col_selected = NULL, col_ID = NULL, cex_all = 0.7,
                          cex_selected = 1, cex_ID = 1, pch_all = 16,
                          pch_selected = 16, mar = NULL) {
  # Initial tests
  if (missing(master)) {
    stop("Argument 'master' is required to produce the plot.")
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

  # Preparing data
  ublocks <- unique(master$data_matrix$Block)
  nblocks <- length(ublocks)

  gvars <- c("Longitude", "Latitude")

  # Colors
  col_pal <- purplow

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

  ## Box to plot
  boxpam <- matrix(terra::ext(master$region), nrow = 2)
  boxpam <- terra::vect(boxpam, crs = terra::crs(master$region))

  ## Plot
  terra::plot(boxpam, col = NA, axes = FALSE, legend = FALSE, mar = mar)
  maps::map(fill = TRUE, col = "gray97", lforce = "n",
            border = "gray80", add = TRUE)
  box(which = "plot")
  points(master$data_matrix[, gvars], pch = pch_all, cex = cex_all,
         col = col_all)

  if (region_border == TRUE) {
    terra::plot(master$region, border = "gray50", add = TRUE)
  }
  if (mask_border == TRUE & !is.null(master$mask)) {
    terra::plot(master$mask, border = "gray50", add = TRUE)
  }

  ## Selected blocks
  if (which == "selected") {
    sel <- which(master$data_matrix$Selected_blocks == 1)
    selected_data <- master$data_matrix[sel, ]
    points(selected_data[sel, gvars], pch = pch_selected, cex = cex_selected,
           col = col_selected)
  }
}
