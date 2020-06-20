#' Selection of survey sites maximizing uniformity in environmental space
#' considering geographic structure
#'
#' @description Selection of sites to be sampled in a survey, with the goal of
#' maximizing uniformity of points in environmental, but considering geographic
#' patterns of data. Similar environments that have a disjoint pattern are selected
#' twice (two survey sites are placed so they consider the biggest geographic
#' clusters).
#'
#' @param master a master_matrix object derived from function
#' \code{\link{master_matrix}} or a master_selection object derived from functions
#' \code{\link{random_selection}}, \code{\link{uniformG_selection}},
#' or \code{\link{uniformE_selection}}.
#' @param variable_1 (character or numeric) name or position of the first
#' variable (X-axis).
#' @param variable_2 (character or numeric) name or position of the second
#' variable (Y-axis).
#' @param select_point (character) How or which point will be selected. Three
#' options are available: "random", "E_centroid", "G_centroid". E_ or G_ centroid
#' indicate that the point(s) closests to the respective centroid will be selected.
#' Default = "E_centroid".
#' @param cluster_method (character) name of the method to be used for detecting
#' clusters. Options are "hierarchical" and "k-means"; default = "hierarchical".
#' See details in \code{\link{find_clusters}}
#' @param sample_for_distance (numeric) sample to be considered when measuring
#' the geographic distances among points in the blocks of environmental points.
#' Default = 250.
#'
#' @return
#' A master_selection object (S3) with an aditional element called
#' selected_sites_EG containing one set of selected sites.
#'
#' @usage
#' EG_selection(master, variable_1, variable_2, select_point = "E_centroid",
#'              cluster_method = "hierarchical", sample_for_distance = 250)
#'
#' @export
#'
#' @examples
#' # Data
#' data("m_matrix", package = "biosurvey")
#'
#' # Making blocks for analysis
#' m_blocks <- make_blocks(m_matrix, variable_1 = "PC1", variable_2 = "PC2",
#'                         n_cols = 10, n_rows = 10, block_type = "equal_area")
#'
#' # Selecting blocks for analysis
#' sel_blocks <- block_sample(m_blocks, variable_1 = "PC1", variable_2 = "PC2",
#'                            expected_blocks = 6, initial_distance = 2,
#'                            increase = 0.1, replicates = 3)
#'
#' # Checking column names
#' colnames(sel_blocks$master_matrix)
#'
#' # Selecting sites uniformly in E and G spaces
#' EG_sel <- EG_selection(master = sel_blocks, variable_1 = "PC1", variable_2 = "PC2",
#'                        select_point = "E_centroid", cluster_method = "hierarchical",
#'                        sample_for_distance = 100)
#'
#' head(EG_sel$selected_sites_EG[[1]])
#' dim(EG_sel$selected_sites_EG[[1]])


EG_selection <- function(master, variable_1, variable_2, select_point = "E_centroid",
                         cluster_method = "hierarchical", sample_for_distance = 250) {

  # Initial tests
  if (missing(master)) {
    stop("Argument 'master' must be defined")
  }
  if (missing(variable_1)) {
    stop("Argument 'variable_1' must be defined.")
  }
  if (missing(variable_2)) {
    stop("Argument 'variable_2' must be defined.")
  }
  if (!select_point[1] %in% c("random", "E_centroid", "G_centroid")) {
    stop("Argument 'select_point' is not valid.")
  }
  if (!cluster_method[1] %in% c("hierarchical", "k-means")) {
    stop("Argument 'cluster_method' is not valid.")
  }
  if (is.null(master$master_matrix$Selected_blocks)) {
    rule <- is.numeric(master$master_matrix$Block)
  } else {
    rule <- master$master_matrix$Selected_blocks == 1
  }

  blocksp <- unique(master$master_matrix[rule, "Block"])
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
  if (nrow(nmodp) > 0) {
    unselp <- point_sample(master$master_matrix[master$master_matrix[, "Block"] %in%
                                                  nmodp[, "Block"], ], variable_1,
                           variable_2, n = 1, select_point = "random",
                           id_column = "Block")
  } else {
    unselp <- matrix(nrow = 0, ncol = ncol(master$master_matrix))
    colnames(unselp) <- colnames(master$master_matrix)
    unselp <- as.data.frame(unselp)
  }

  # unimodal
  if (nrow(unimp) > 0) {
    ueselp <- point_sample(master$master_matrix[master$master_matrix[, "Block"] %in%
                                                  unimp[, "Block"], ], variable_1,
                           variable_2, n = 1, select_point = select_point,
                           id_column = "Block")
  } else {
    ueselp <- matrix(nrow = 1, ncol = ncol(master$master_matrix))
    colnames(ueselp) <- colnames(master$master_matrix)
    ueselp <- as.data.frame(na.omit(ueselp))
  }

  # multimodal
  if (nrow(mmodp) > 0) {
    distsp <- distsp[as.character(mmodp$Block)]
    meselp <- point_sample_cluster(master$master_matrix[master$master_matrix[, "Block"] %in%
                                                          mmodp[, "Block"], ],
                                   variable_1, variable_2, distance_list = distsp,
                                   n = 1, cluster_method = cluster_method,
                                   select_point = select_point, id_column = "Block")
  } else {
    meselp <- matrix(nrow = 0, ncol = ncol(master$master_matrix))
    colnames(meselp) <- colnames(master$master_matrix)
    meselp <- as.data.frame(meselp)
  }

  # combining
  master$selected_sites_EG <- list(selection_1 = rbind(unselp, ueselp, meselp))

  return(structure(master, class = "master_selection"))
}
