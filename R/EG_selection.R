#' Selection of survey sites maximizing uniformity in both environmental and
#' geographic space
#'
#' @description Selection of sites to be sampled in a survey, with the goal of
#' maximizing uniformity of points in environmental and geographic space.
#'
#' @param master a master_matrix object derived from function
#' \code{\link{master_matrix}} or a master_selection object derived from functions
#' \code{\link{random_selection}}, \code{\link{uniformG_selection}},
#' or \code{\link{E_selection}}.
#' @param variable_1 (character or numeric) name or position of the first
#' variable (X-axis).
#' @param variable_2 (character or numeric) name or position of the second
#' variable (Y-axis).
#' @param select_point (character) Default = "E_centroid".
#' @param cluster_method (character) Default = "hierarchical".
#' @param sample_for_distance (numeric) Default = 250.
#'
#' @return
#' A master_selection object (S3) with an aditional element called
#' selected_sites_EG containing one set of selected sites.
#'
#' @usage
#' EG_selection(master, variable_1, varaible_2, select_point = "E_centroid",
#'              cluster_method = "hierarchical", sample_for_distance = 250)
#'
#' @export


EG_selection <- function(master, variable_1, varaible_2, select_point = "E_centroid",
                         cluster_method = "hierarchical", sample_for_distance = 250) {

  # Initial tests
  if (missing(master)) {
    stop("Argument 'master' must be defined")
  }
  if (missing(variable_1)) {
    stop("Argument 'variable_1' must be defined.")
  }
  if (missing(varaible_2)) {
    stop("Argument 'varaible_2' must be defined.")
  }
  if (!select_point[1] %in% c("random", "E_centroid", "G_centroid")) {
    stop("Argument 'select_point' is not valid.")
  }
  if (!cluster_method[1] %in% c("hierarchical", "k-means")) {
    stop("Argument 'cluster_method' is not valid.")
  }

  blocksp <- unique(master$master_matrix[, "Block"])
  g_cols <- c("Longitude", "Latitude")

  distsp <- lapply(blocksp, function(x) {
    block_data <- master$master_matrix[master$master_matrix[, "Block"] == x,
                                       g_cols]
    tpoints <- nrow(block_data)

    if (tpoints > sample_for_distance) {
      block_data <- block_data[sample(tpoints, sample_for_distance), ]
    }

    dsnna <- na.omit(c(raster::pointDistance(block_data, lonlat = TRUE)))
    dsnna[dsnna != 0]
  })
  names(distsp) <- blocksp

  dpp <- unimodal_test(distsp)

  unimp <- dpp[which(dpp$p_alue > 0.05), ]
  mmodp <- dpp[which(dpp$p_alue <= 0.05), ]
  nmodp <- dpp[which(is.na(dpp$p_alue)), ]

  # no mode (very few points)
  unselp <- point_sample(master$master_matrix[master$master_matrix[, "Block"] %in%
                                                nmodp[, "Block"], ], variable_1,
                         varaible_2, n = 1, select_point = "random",
                         id_column = "Block")

  # unimodal
  ueselp <- point_sample(master$master_matrix[master$master_matrix[, "Block"] %in%
                                                unimp[, "Block"], ], variable_1,
                         varaible_2, n = 1, select_point = select_point,
                         id_column = "Block")

  # multimodal
  meselp <- point_sample_cluster(master$master_matrix[master$master_matrix[, "Block"] %in%
                                                        mmodp[, "Block"], ],
                                 variable_1, varaible_2, distance_list = distsp,
                                 n = 1, cluster_method = cluster_method,
                                 select_point = select_point, id_column = "Block")

  # combining
  master$selected_sites_EG <- list(selection_1 = rbind(ueselp, unselp, meselp))

  return(structure(master, class = "master_selection"))
}
