#' Creates a block-like regionalization of environmental space
#'
#' @description Divides a two-dimensional cloud of points in blocks according to
#' a user-defined number of rows and columns. This is applied to the element
#' master_matrix and, if not NULL, to preselected_sites.
#'
#' @param master_matrix object derived from function \code{\link{prepare_master_matrix}}.
#' @param variable_1 (character or numeric) name or position of the first
#' variable (X axis) to be used to create blocks.
#' @param variable_2 (character or numeric) name or position of the second
#' variable (Y axis) to be used to create blocks (must be different from the
#' first one).
#' @param n_cols (numeric) number of columns of a grid used to creates blocks and
#' split the bi-dimensional space.
#' @param n_rows (numeric) number of rows of a grid used to creates blocks and
#' split the bi-dimensional space. If NULL, the default, \code{n_cols = n_rows}.
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
#' containing block identifiers. If the element preselected_sites is not NULL in
#' master_matrix blocks are also assigned to this sites.
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
#' unique(m_blocks$data_matrix$Block)


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
  coln <- colnames(master_matrix$data_matrix)
  if (!variable_1 %in% coln) {
    stop(variable_1, " is not one o the columns in 'master_matrix$data_matrix'.")
  }
  if (!variable_2 %in% coln) {
    stop(variable_2, " is not one o the columns in 'master_matrix$data_matrix'.")
  }

  # Preparing data
  data <- master_matrix$data_matrix
  id <- paste(data[, 1], data[, 2])

  if (!is.null(master_matrix$preselected_sites)) {
    predata <- master_matrix$preselected_sites
    idpre <- paste(predata[, 2], predata[, 3])
  }

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
    all_cls <- assign_blocks(data, variable_1, variable_2, n_cols, n_rows, xlb,
                             ylb, block_type = "equal_area")

    # Assigning blocks to user predefined sites
    if (!is.null(master_matrix$preselected_sites)) {
      prese <- assign_blocks(predata, variable_1, variable_2, n_cols, n_rows, xlb,
                             ylb, block_type = "equal_area")
    }
  } else {
    # Detecting ranges and intervals
    xlb <- seq(0, 1, (1 / n_cols))
    xlb[length(xlb)] <- 1

    # Assigning block numbers
    all_cls <- assign_blocks(data, variable_1, variable_2, n_cols, n_rows, xlb,
                             block_type = "equal_points")

    # Assigning blocks to user predefined sites
    if (!is.null(master_matrix$preselected_sites)) {
      prese <- assign_blocks(predata, variable_1, variable_2, n_cols, n_rows, xlb,
                             block_type = "equal_points")
    }
  }

  # Returning results
  all_cls <- all_cls[match(id, paste(all_cls[, 1], all_cls[, 2])), ] # matches data back in order
  master_matrix$data_matrix <- all_cls

  if (!is.null(master_matrix$preselected_sites)) {
    prese <- prese[match(idpre, paste(prese[, 2], prese[, 3])), ] # matches preselected data back in order
    master_matrix$preselected_sites <- prese
  }

  return(structure(master_matrix, class = "master_matrix"))
}



#' Helper to assign block numbers to data according to variables and limits
#'
#' @param data a matrix or a data frame that contains at least four columns:
#' "Longitude" and "Latitude" to represent geographic position, and two other
#' columns to represent the variables of the 2D environmental space.
#' @param variable_1 (character or numeric) name or position of the first
#' variable (X axis) to be used to create blocks.
#' @param variable_2 (character or numeric) name or position of the second
#' variable (Y axis) to be used to create blocks (must be different from the
#' first one).
#' @param n_cols (numeric) number of columns of a grid used to creates blocks and
#' split the bi-dimensional space.
#' @param n_rows (numeric) number of rows of a grid used to creates blocks and
#' split the bi-dimensional space. If NULL, the default, \code{n_cols = n_rows}.
#' @param xlb (numeric) Vector of values of extremes for all blocks considering
#' \code{variable_1}.
#' @param ylb (numeric) Vector of values of extremes for all blocks considering
#' \code{variable_2}. Needed when \code{block_type} = "equal area". Default = NULL.
#' @param block_type (character) type of blocks to be use for dividing
#' the bi-dimensional space. Two options are available: "equal_area" and
#' "equal_points". Default = "equal_area".
#'
#' @return
#' Original element defined in \code{data} plus a new column named "Block"
#' defining the block that correspond to each of the points represented in rows.
#'
#' @export
#'
#' @usage
#' assign_blocks(data, variable_1, variable_2, n_cols, n_rows = NULL,
#'               xlb, ylb = NULL, block_type = "equal_area")
#' @examples
#' # data
#' dat <- matrix(runif(800), ncol = 4)
#' xlims <- quantile(dat[, 3])
#' ylims <- quantile(dat[, 4])
#'
#' # assigning blocks
#' datb <- assign_blocks(dat, variable_1 = 3, variable_2 = 4, n_cols = 10,
#'                       xlb = xlims, ylb = ylims, block_type = "equal_area")

assign_blocks <- function(data, variable_1, variable_2, n_cols,
                          n_rows = NULL, xlb, ylb = NULL,
                          block_type = "equal_area") {
  # Initial tests
  if (missing(data)) {
    stop("Argument 'data' needs to be defined.")
  }
  if (missing(variable_1)) {
    stop("Argument 'variable_1' needs to be defined.")
  }
  if (missing(variable_2)) {
    stop("Argument 'variable_2' needs to be defined.")
  }
  if (missing(n_cols)) {
    stop("Argument 'n_cols' needs to be defined.")
  }
  if (is.null(n_rows)) {
    n_rows <- n_cols
  }
  if (missing(xlb)) {
    stop("Argument 'xlb' needs to be defined.")
  }
  if (is.null(ylb) & block_type == "equal_area") {
    stop("Argument 'ylb' needs to be defined.")
  }
  if (!block_type[1] %in% c("equal_area", "equal_points")) {
    stop("Argument 'block_type' is not valid.")
  }


  # assigning blocks
  if (block_type == "equal_area") {
    ## blocks of equal area
    all_cls <- lapply(1:(length(xlb) - 1), function(x) {
      ## x axis
      if(x == 1){
        x1 <- data[, variable_1] >= xlb[x]
      } else {
        x1 <- data[, variable_1] > xlb[x]
      }
      xid <- which(x1 & data[, variable_1] <= xlb[(x + 1)])
      if (length(xid) > 0) {
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
      } else {
        return(data[0, ])
      }
    })

    # Finishing assigning
    all_cls <- do.call(rbind, all_cls)
    colnames(all_cls)[ncol(all_cls)] <- "Block"
    all_cls <- all_cls[order(all_cls[, "Block"]), ]
  } else {
    ## blocks with equal number of points
    all_cls <- lapply(1:(length(xlb) - 1), function(x) {
      ## x axis
      q1 <- quantile(data[, variable_1], xlb[x])
      x1 <- ifelse(x == 1, q1, (q1 + 0.000000000000001))
      x2 <- quantile(data[, variable_1], xlb[(x + 1)])
      xid <- which(data[, variable_1] >= x1 & data[, variable_1] <= x2)
      if (length(xid) > 0) {
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
      } else {
        return(data[0, ])
      }

    })

    # Finishing assigning
    all_cls <- do.call(rbind, all_cls)
    colnames(all_cls)[ncol(all_cls)] <- "Block"
  }
  return(all_cls)
}
