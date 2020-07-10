#' Creates a block-like regionalization of environmental space
#'
#' @description Divides a two-dimensional cloud of points in blocks according to
#' a user-defined number of rows and columns.
#'
#' @param master_matrix object derived from function \code{\link{master_matrix}}.
#' @param variable_1 (character or numeric) name or position of the first
#' variable (X axis) to be used to create blocks.
#' @param variable_2 (character or numeric) name or position of the second
#' variable (Y axis) to be used to create blocks (must be different from the
#' first one).
#' @param n_cols (numeric) number of columns of a grid used to creates blocks and
#' split the bi-dimensional space.
#' @param n_rows (numeric) number of rows of a grid used to creates blocks and
#' split the bi-dimensional space. If not defined, \code{n_cols = n_rows}.
#' Default = NULL.
#' @param block_type (character) type of blocks to be use for dividing
#' the bi-dimensional space. Two options are available: "equal_area" and
#' "equal_points". Default = "equal_area".
#'
#' @details
#' For block_type, option "equal_area" generates blocks of the same size. The other
#' option ("equal_points"), generates blocks containing the same number of points,
#' which generally results in blocks of different sizes.
#'
#' @return
#' An S3 object of class master_matrix, containing the same elements found in a
#' master_matrix object, with an additional column on the master_matrix data.frame
#' containing block identifiers.
#'
#' @usage
#' make_blocks(master_matrix, variable_1, variable_2, n_cols, n_rows = NULL,
#'             block_type = "equal_area")
#'
#' @export
#' @importFrom stats quantile
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#'
#' # Creating blocks
#' m_blocks <- make_blocks(m_matrix, variable_1 = "PC1",
#'                         variable_2 = "PC2", n_cols = 10, n_rows = 10,
#'                         block_type = "equal_area")
#' unique(m_blocks$master_matrix$Block)


make_blocks <- function(master_matrix, variable_1, variable_2, n_cols,
                        n_rows = NULL, block_type = "equal_area") {
  # Initial tests
  if (missing(master_matrix)) {
    stop("Argument 'master_matrix' needs to be defined.")
  }
  if (missing(variable_1) | missing(variable_2)) {
    stop("Argument 'variable_1' and 'variable_2' needs to be defined.")
  }
  if (missing(n_cols)) {
    stop("Argument 'n_cols' needs to be defined.")
  }
  if (is.null(n_rows)) {
    n_rows <- n_cols
  }
  if (!block_type[1] %in% c("equal_area", "equal_points")) {
    stop("Argument 'block_type' is not valid.")
  }

  # Preparing data
  data <- master_matrix$master_matrix
  id <- paste(data[, 1], data[, 2])

  # Block partition
  if (block_type[1] == "equal_area") {
    # Detecting ranges and intervals
    xrange <- range(data[, variable_1])
    xinter <- diff(xrange) / n_cols
    yrange <- range(data[, variable_2])
    yinter <- diff(yrange) / n_rows

    xlb <- seq(xrange[1], xrange[2], xinter)
    xlb[length(xlb)] <- xrange[2]
    ylb <- seq(yrange[1], yrange[2], yinter)
    ylb[length(ylb)] <- yrange[2]

    # Assigning block numbers
    all_cls <- lapply(1:(length(xlb) - 1), function(x) {
      ## x axis
      if(x == 1){
        x1 <- data[, variable_1] >= xlb[x]
      } else {
        x1 <- data[, variable_1] > xlb[x]
      }
      xid <- which(x1 & data[, variable_1] <= xlb[(x + 1)])
      pd <- data[xid, ]
      pd <- cbind(pd, NA)
      if (nrow(pd) > 0) {
        ## y axis
        for (y in 1:(length(ylb) - 1)) {
          if(y == 1) {
            y1 <- pd[, variable_2] >= ylb[y]
          } else {
            y1 <- pd[, variable_2] > ylb[y]
          }
          yid <- which(y1 & pd[, variable_2] <= ylb[(y + 1)])
          nb <- ifelse(x == 1, y, (x * length(ylb)) + y)
          pd[yid, ncol(pd)] <- rep(nb, length(yid))
        }
      }
      return(pd)
    })

    # Finishing assigning
    all_cls <- do.call(rbind, all_cls)
    colnames(all_cls)[ncol(all_cls)] <- "Block"
    all_cls <- all_cls[order(all_cls[, "Block"]), ]
    ub <- unique(all_cls[, "Block"])
    blks <- lapply(1:length(ub), function(x) {
      rep(x, sum(all_cls[, "Block"] == ub[x]))
    })
    all_cls[, "Block"] <- unlist(blks)

  } else {
    # Detecting ranges and intervals
    xlb <- seq(0, 1, (1 / n_cols))
    xlb[length(xlb)] <- 1

    # Assigning block numbers
    all_cls <- lapply(1:(length(xlb) - 1), function(x) {
      ## x axis
      q1 <- quantile(data[, variable_1], xlb[x])
      x1 <- ifelse(x == 1, q1, (q1 + 0.000000000000001))
      x2 <- quantile(data[, variable_1], xlb[(x + 1)])
      xid <- which(data[, variable_1] >= x1 & data[, variable_1] <= x2)
      pd <- data[xid, ]
      pd <- cbind(pd, NA)

      if (n_cols != n_rows) {
        ylb <-  seq(0, 1, round(1 / n_rows, 5))
        ylb[length(ylb)] <- 1
      } else {
        ylb <- xlb
      }
      ## y axis
      for (y in 1:(length(ylb) - 1)) {
        qy1 <- quantile(pd[, variable_2], ylb[y])
        y1 <- ifelse(y == 1, qy1, (qy1 + 0.000000000000001))
        y2 <- quantile(pd[, variable_2], ylb[(y + 1)])
        yid <- which(pd[, variable_2] >= y1 & pd[, variable_2] <= y2)
        nb <- ifelse(x == 1, y, ((x - 1) * (length(ylb) - 1)) + y)
        pd[yid, ncol(pd)] <- rep(nb, length(yid))
      }
      return(pd)
    })

    # Finishing assigning
    all_cls <- do.call(rbind, all_cls)
    colnames(all_cls)[ncol(all_cls)] <- "Block"
  }

  # Returning results
  all_cls <- all_cls[match(id, paste(all_cls[, 1], all_cls[, 2])), ]
  master_matrix$master_matrix <- all_cls
  return(structure(master_matrix, class = "master_matrix"))
}
