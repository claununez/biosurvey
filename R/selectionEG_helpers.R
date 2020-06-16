#'
#'
#' @description
#'
#' @param data a matrix or a data frame that contains at least two columns.
#' @param variable_1 (character or numeric) name or position of the first
#' variable (X axis).
#' @param variable_2 (character or numeric) name or position of the second
#' variable (Y axis). Must be different from the first one.
#' @param n (numeric) number of points that are close to the centroid to be
#' detected. Default = 1.
#' @param select_point (character) Three options are available: "random",
#' "E_centroid", "G_centroid". Default = "E_centroid".
#' @param id_column (character or numeric) name or numeric index of the column
#' in \code{data} containing identifiers of one or distint sets of points.
#' If, NULL, the default, only one set is assumed.
#'
#' @return
#'
#' @usage
#' point_sample(data, variable_1, varaible_2, n = 1, select_point = "E_centroid",
#'              id_column = NULL)
#'
#' @export
#'

point_sample <- function(data, variable_1, varaible_2, n = 1,
                         select_point = "E_centroid", id_column = NULL) {
  # initial tests
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (missing(variable_1)) {
    stop("Argument 'variable_1' must be defined.")
  }
  if (missing(varaible_2)) {
    stop("Argument 'varaible_2' must be defined.")
  }
  if (!select_point[1] %in% c("random", "E_centroid", "G_centroid")) {
    stop("Argument 'select_point' is not valid, options are:\n'random', 'E_centroid', 'G_centroid'")
  }

  e_cols <- c(variable_1, variable_2)
  g_cols <- c("Longitude", "Latitude")
  bda <- data[, id_column]
  bs <- unique(bda)

  if (select_point[1] == "random") {
    samp <- lapply(bs, function(x) {
      bd <- data[bda == x, ]
      bd <- bd[sample(1:nrow(bd), n), ]
      unique(bd)
    })
    bsam <- do.call(rbind, samp)
    colnames(bsam) <- colnames(data)
  }

  if (select_point[1] == "E_centroid") {
    bsam <- closest_to_centroid(data, e_cols[1], e_cols[2], space = "E",
                                id_column = id_column)
  }

  if (select_point[1] == "G_centroid") {
    bsam <- closest_to_centroid(data, g_cols[1], g_cols[2], space = "G",
                                id_column = id_column)
  }

  return(bsam)
}





#'
#'
#' @description
#'
#' @param density
#'
#' @return
#'
#' @usage
#' find_modes(density)
#'
#' @export
#'

find_modes <- function(density) {

  # initial tests
  if (missing(density)) {
    stop("Argument 'density' must be defined.")
  }

  density_y <- density$y
  modes <- NULL

  for (i in 2:(length(density_y) - 1)) {
    if ((density_y[i] > density_y[i - 1]) & (density_y[i] > density_y[i + 1])) {
      modes <- c(modes, i)
    }
  }

  if ( length(modes) == 0 ) {
    modes <- "This is a monotonic distribution."
    return(modes)
  } else {
    return(data.frame(mode = density$x[modes], density = density_y[modes]))
  }
}





#'
#'
#' @description
#'
#' @param data a matrix or a data frame that contains at least two columns.
#' @param x_column (character) the name of the X-axis.
#' @param y_column (character) the name of the Y-axis.
#' @param space (character) space in which the thinning will be performed. There
#' are two options available: "G", if it will be in the geographyc space, and
#' "E", if it will be on the environmental space.
#' @param cluster_method (character) There are two options available:
#' "hierarchical" and "k-means". Default = "hierarchical".
#' @param n_kmeans Default = NULL.
#' @param split_distance Default = NULL.
#'
#' @return
#'
#' @usage
#' find_clusters(data, x_column, y_column, space,
#'               cluster_method = "hierarchical", n_kmeans = NULL,
#'               split_distance = NULL)
#'
#' @export
#' @importFrom stats hclust cutree kmeans dist
#' @importFrom raster::pointDistance
#'

find_clusters <- function(data, x_column, y_column, space, cluster_method = "hierarchical",
                          n_kmeans = NULL, split_distance = NULL) {

  # initial tests
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (missing(x_column)) {
    stop("Argument 'x_column' must be defined.")
  }
  if (missing(y_column)) {
    stop("Argument 'y_column' must be defined.")
  }
  if (missing(space)) {
    stop("Argument 'space' is not defined.")
  }

  if (cluster_method %in% c("hierarchical", "k-means")) {
    if (cluster_method[1] == "hierarchical") {
      is.null(split_distance) {
        stop("Argument 'split_distance' must be defined if 'cluster_method' = 'hierarchical'.")
      }

      if (space == "E") {
        cluster <- stats::hclust(dist(data[, c(x_column, y_column)]),
                                 method = "complete")
      } else {
        cluster <- stats::hclust(as.dist(raster::pointDistance(data[, c(x_column, y_column)],
                                                               lonlat = T)),
                          method = "complete")
      }

      cluster_vector <- stats::cutree(cluster, h = split_distance)

    } else {
      is.null(n_kmeans) {
        stop("Argument 'n_kmeans' must be defined if 'cluster_method' = 'k-means'.")
      }

      set.seed(1)
      cluster_vector <- stats::kmeans(as.matrix(data[, c(x_column, y_column)]),
                               n_kmeans)$cluster
    }
  } else {
    stop("Argument 'cluster_method' is not valid.")
  }

  data <- data.frame(data, clusters = cluster_vector)
  return(data)
}





#'
#'
#' @description
#'
#' @param data a matrix or a data frame that contains at least two columns.
#' @param variable_1 (character or numeric) name or position of the first
#' variable (X axis).
#' @param variable_2 (character or numeric) name or position of the second
#' variable (Y axis). Must be different from the first one.
#' @param distance_list
#' @param n (numeric) number of points that are close to the centroid to be
#' detected. Default = 1.
#' @param cluster_method (character) There are two options available:
#' "hierarchical" and "k-means". Default = "hierarchical".
#' @param select_point (character) Three options are available: "random",
#' "E_centroid", "G_centroid". Default = "E_centroid".
#' @param id_column (character or numeric) name or numeric index of the column
#' in \code{data} containing identifiers of one or distint sets of points.
#' If, NULL, the default, only one set is assumed.
#'
#' @return
#'
#' @usage
#' point_sample_cluster(data, variable_1, varaible_2, distance_list,
#'                      n = 1, cluster_method = "hierarchical",
#'                      select_point = "E_centroid", id_column = NULL)
#'
#' @export
#'

point_sample_cluster <- function(data, variable_1, varaible_2, distance_list,
                                 n = 1, cluster_method = "hierarchical",
                                 select_point = "E_centroid", id_column = NULL) {
  # initial tests
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (missing(variable_1)) {
    stop("Argument 'variable_1' must be defined.")
  }
  if (missing(varaible_2)) {
    stop("Argument 'varaible_2' must be defined.")
  }
  if (!select_point[1] %in% c("random", "E_centroid", "G_centroid")) {
    stop("Argument 'select_point' is not valid, options are:\n'random', 'E_centroid', 'G_centroid'")
  }
  if (missing(distance_list)) {
    stop("Argument 'distance_list' must be defined")
  }
  if (!cluster_method[1] %in% c("hierarchical", "k-means")) {
    stop("Argument 'cluster_method' is not valid.")
  }

  bda <- data[, id_column]
  bs <- unique(bda)

  mgsel <- lapply(bs, function(x) {
    md <- find_modes(density = density(distance_list[[as.character(x)]]))

    if (nrow(md) > 1) {
      dens <- md$density
      mdss <- md[order(dens), ]
      dd <- abs(diff(mdss[(length(dens) - 1):length(dens), 1]))

      clush <- find_clusters(data = data[bda == x, ], variable_1, variable_2,
                             space = "G", cluster_method = cluster_method,
                             split_distance = dd)

      sel <- as.numeric(names(sort(table(clush$clusters), decreasing = T)[1:2]))

      bse <- point_sample(data = clush[clush$clusters %in% sel, ],
                          variable_1, varaible_2, n = n,
                          select_point = select_point, id_column = id_column)
      bse$clusters <- NULL
    } else {
      bse <- point_sample(data = data[bda == x, ], variable_1, varaible_2, n = n,
                          select_point = select_point, id_column = id_column)
    }
    return(bse)
  })
  mgsel <- do.call(rbind, mgsel)

  return(mgsel)
}





#'
#'
#' @description
#'
#' @param distance_list
#'
#' @return
#'
#' @usage
#' unimodal_test(distance_list)
#'
#' @export
#' @importFrom diptest dip.test
#'
unimodal_test <- function(distance_list) {

  # initial tests
  if (missing(distance_list)) {
    stop("Argument 'distance_list' must be defined.")
  }

  bs <- names(distance_list)

  dss <- lapply(bs, function(x) {
    ds <- distance_list[[x]]

    if (length(ds) <= 2) {
      return(data.frame(Block = x, D = NA, p_alue = NA))
    } else {
      dp <- diptest::dip.test(ds, simulate.p.value = TRUE, B = 1000)
      return(data.frame(Block = x, D = dp$statistic, p_alue = dp$p.value))
    }
  })

  return(do.call(rbind, dss))
}
