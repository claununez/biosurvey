#' Helps in thinning points either in geographic or environmental space
#'
#' @description Point thinning based on user-defined distances in geographic or
#' environmental space.
#'
#' @param data a matrix or a data frame that contains at least two columns.
#' @param x_column (character) the name of the X-axis.
#' @param y_column (character) the name of the Y-axis.
#' @param thinning_distance (numeric) distance for thinning. Units must be
#' selected according to the space, km for geographic and euclidean distances for
#' environmental space.
#' @param space (character) space in which the thinning will be performed. There
#' are two options available: "G", if it will be in the geographyc space, and
#' "E", if it will be on the environmental space.
#' @param max_n_samples (numeric) maximun number of samples to chose with most
#' points included. Default = 1.
#' @param replicates (numeric) number of thinning replicates. Default = 10.
#' @param set_seed (numeric) integer value to specify a initial seed. Default = 1.
#'
#' @return
#' A list with one or more elements, depending on \code{max_n_samples}. Each
#' element is a data.frame containing points retained after thinning. All elements
#' are different in at least one of the selected points.
#'
#' @usage
#' point_thinning(data, x_column, y_column, thinning_distance, space,
#'                max_n_samples = 1, replicates = 10, set_seed = 1)
#'
#' @export
#' @importFrom sp coordinates
#' @importFrom spatstat ppp closepairs
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#' data1 <- m_matrix$master_matrix
#'
#' # Thinnig the points
#' thin <- point_thinning(data1, x_column = "Longitude", y_column = "Latitude",
#'                        thinning_distance = 200, space = "G",
#'                        max_n_samples = 1, replicates = 5, set_seed = 1)

point_thinning <- function(data, x_column, y_column, thinning_distance, space,
                           max_n_samples = 1, replicates = 10, set_seed = 1) {
  # Initial tests
  if (missing(data)) {
    stop("Argument 'data' is not defined.")
  }
  if (missing(x_column)) {
    stop("Argument 'x_column' is not defined.")
  }
  if (missing(y_column)) {
    stop("Argument 'y_column' is not defined.")
  }
  if (missing(thinning_distance)) {
    stop("Argument 'thinning_distance' is not defined.")
  }
  if (missing(space)) {
    stop("Argument 'space' is not defined.")
  }

  # Initial preparation
  data <- data[!is.na(data[, x_column]) & !is.na(data[, y_column]), ]
  data <- data[!duplicated(paste(data[, x_column], data[, y_column], sep = "_")), ]

  cls <- class(data)[1]
  if (cls != "data.frame") {
    if (cls == "matrix") {
      data <- as.data.frame(data)
    } else {
      stop("'data' must be of class matrix or data.frame.")
    }
  }

  # Preprocessing if space = G
  if (space == "G") {
    data_sp <- wgs84_2aed_laea(data, x_column, y_column, which = "ED")
    xy <- sp::coordinates(data_sp)
    data$xaed <- xy[, 1]
    data$yaed <- xy[, 2]
    xyo <- c(x_column, y_column)
    x_column <- "xaed"
    y_column <- "yaed"
    thinning_distance <- thinning_distance * 1000
  }

  # Thinning points
  lthins <- lapply(1:replicates, function(x) {
    set.seed(set_seed + x - 1)

    data1 <- data[sample(nrow(data)), ]
    X <- spatstat::ppp(data1[, x_column], data1[, y_column],
                       range(data1[, x_column]), range(data1[, y_column]))

    close_index <- spatstat::closepairs(X, rmax = thinning_distance,
                                        what = "indices", twice = FALSE)$j

    return(data1[!1:nrow(data1) %in% close_index, !colnames(data1) %in%
                   c("xaed", "yaed")])
  })

  if (space == "G") {
    x_column <- xyo[1]; y_column <- xyo[2]
  }

  # Getting needed samples form the most numerous ones
  lns <- sapply(lthins, nrow)
  lthins <- lthins[which(lns == max(lns))]

  cd <- sapply(lthins, function(x) {
    y <- x[order(x[, x_column], x[, y_column]), ]
    paste0(paste0(y[, x_column], y[, y_column]), collapse = "_")
  })

  lthins <- lthins[which(!duplicated(cd))]
  nsel <- ifelse(length(lthins) < max_n_samples, length(lthins), max_n_samples)

  # Returning results
  if (nsel == 1) {
    return(list(lthins[[1]]))
  } else {
    return(lthins[1:nsel])
  }
}





#' Detection of the closest points to the centroid of a cloud of points
#'
#' @param data a matrix or a data frame that contains at least two columns.
#' @param x_column (character) the name of the X-axis.
#' @param y_column (character) the name of the Y-axis.
#' @param space (character) space in which the thinning will be performed. There
#' are two options available: "G", if it will be in the geographyc space, and
#' "E", if it will be on the environmental space.
#' @param n (numeric) number of points that are close to the centroid to be
#' detected. Default = 1.
#' @param id_column (character or numeric) name or numeric index of the column
#' in \code{data} containing identifiers of one or distint sets of points.
#' If, NULL, the default, only one set is assumed.
#'
#' @return
#' A data.frame containing \code{n} rows corresponding to the point or points that
#' are the closest to the cenroid of all other points of reference.
#'
#' @usage
#' closest_to_centroid(data, x_column, y_column, space, n = 1, id_column = NULL)
#'
#' @export
#' @importFrom raster pointDistance
#' @importFrom stats mahalanobis qchisq cov dist
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#' data1 <- m_matrix$master_matrix
#'
#' # Finding the closest point to the centroid
#' centroid <- closest_to_centroid(data1, x_column = "Longitude",
#'                                 y_column = "Latitude", space = "G",
#'                                 n = 1, id_column = NULL)

closest_to_centroid <- function(data, x_column, y_column, space, n = 1,
                                id_column = NULL) {
  # Initial tests
  if (missing(data)) {
    stop("Argument 'data' is not defined.")
  }
  if (missing(x_column)) {
    stop("Argument 'x_column' is not defined.")
  }
  if (missing(y_column)) {
    stop("Argument 'y_column' is not defined.")
  }
  if (missing(space)) {
    stop("Argument 'space' is not defined.")
  }

  # Detection of sets of points
  if(!is.null(id_column)) {
    bda <- data[, id_column]
    bs <- sort(unique(bda))
  } else {
    data <- cbind(data, id_column = 1)
    id_column <- "id_column"
    bda <- data[, id_column]
    bs <- sort(unique(bda))
  }

  # Detection of points closest to centroid
  ucent <- lapply(bs, function(x) {
    gblock <- data[bda == x, ]
    if (nrow(gblock) > 2) {
      cent <- apply(gblock[, c(x_column, y_column)], 2, mean)
      level <- 0.01
      if (nrow(gblock) > 20) {
        covm <- stats::cov(gblock[, c(x_column, y_column)])
        ndim <- length(cent); sigma_i <- solve(covm) / stats::qchisq(level, df = ndim)
        stds <- 1 / sqrt(eigen(sigma_i)$values)
        hl <- cent + stds; ll <- cent - stds
        c1 <- gblock[, x_column] >= ll[1] & gblock[, x_column] <= hl[1] &
          gblock[, y_column] >= ll[2] & gblock[, y_column] <= hl[2]
        con <- sum(c1)

        if (con <= 2) {
          while (con == 0) {
            if (level > 0.98) {
              if (space == "G") {
                ds <- raster::pointDistance(cent, gblock[, c(x_column, y_column)],
                                            lonlat = TRUE)
              } else {
                ds <- stats::mahalanobis(x = gblock[, c(x_column, y_column)],
                                         center = cent, cov = covm, tol = 0.0000009)
              }
              break()
            }
            level <- level + 0.01
            sigma_i <- solve(covm) / stats::qchisq(level, df = ndim)
            stds <- 1 / sqrt(eigen(sigma_i)$values)
            hl <- cent + stds
            ll <- cent - stds
            c1 <- gblock[, x_column] >= ll[1] & gblock[, x_column] <= hl[1] &
              gblock[, y_column] >= ll[2] & gblock[, y_column] <= hl[2]
            con <- sum(c1)
            if (con > 0) {
              break()
            }
          }
        }
      } else {
        c1 <- rep(TRUE, nrow(gblock))
      }
      if (level > 0.98) {
        return(gblock[which(ds == sort(ds)[1:n])[1:n], ])
      }
      if (space == "G") {
        ds <- raster::pointDistance(cent, gblock[c1, c(x_column, y_column)],
                                    lonlat = TRUE)
        return(gblock[which(ds %in% sort(ds)[1:n])[1:n], ])
      } else {
        ds <- as.matrix(stats::dist(rbind(cent, gblock[c1, c(x_column,
                                                             y_column)])))[-1, 1]
        return(gblock[which(ds %in% sort(ds)[1:n])[1:n], ])
      }
    } else {
      return(gblock[1, ])
    }
  })
  return(do.call(rbind, ucent))
}
