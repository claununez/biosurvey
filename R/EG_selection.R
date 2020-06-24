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
#' @param n_blocks (numeric) number of blocks to be selected from all existent
#' blocks in \code{master$master_matrix}.
#' @param initial_distance (numeric) euclidean distance to be used for a first
#' process of thinning and detection of remaining points.
#' @param increase (numeric) value to be added to \code{initial_distance} until
#' reaching the number of \code{n_blocks}.
#' @param replicates (numeric) number of thinning replicates performed to select
#' blocks uniformly. Default = 10.
#' @param max_n_samplings (numeric) maximum number of samples to be chosen after
#' performing all thinning \code{replicates}. Default = 1.
#' @param select_point (character) How or which point will be selected. Three
#' options are available: "random", "E_centroid", "G_centroid". E_ or G_ centroid
#' indicate that the point(s) closets to the respective centroid will be selected.
#' Default = "E_centroid".
#' @param cluster_method (character) name of the method to be used for detecting
#' clusters. Options are "hierarchical" and "k-means"; default = "hierarchical".
#' See details in \code{\link{find_clusters}}
#' @param sample_for_distance (numeric) sample to be considered when measuring
#' the geographic distances among points in the blocks of environmental points.
#' Default = 250.
#' @param set_seed (numeric) integer value to specify a initial seed. Default = 1.
#'
#' @return
#' A master_selection object (S3) with an additional element called
#' selected_sites_EG containing one or more sets of selected sites depending on
#' \code{max_n_samplings}.
#'
#' @usage
#' EG_selection(master, variable_1, variable_2, n_blocks, initial_distance,
#'              increase, replicates = 10, max_n_samplings = 1,
#'              select_point = "E_centroid", cluster_method = "hierarchical",
#'              sample_for_distance = 250, set_seed = 1)
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
#' # Checking column names
#' colnames(m_blocks$master_matrix)
#'
#' # Selecting sites uniformly in E and G spaces
#' EG_sel <- EG_selection(master = m_blocks, variable_1 = "PC1", variable_2 = "PC2",
#'                        n_blocks = 10, initial_distance = 1.5, increase = 0.1,
#'                        replicates = 1, max_n_samplings = 1,
#'                        select_point = "E_centroid",
#'                        cluster_method = "hierarchical",
#'                        sample_for_distance = 100)
#'
#' head(EG_sel$selected_sites_EG[[1]])
#' dim(EG_sel$selected_sites_EG[[1]])


EG_selection <- function(master, variable_1, variable_2, n_blocks, initial_distance,
                         increase, replicates = 10, max_n_samplings = 1,
                         select_point = "E_centroid", cluster_method = "hierarchical",
                         sample_for_distance = 250, set_seed = 1) {

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
  if (missing(n_blocks)) {
    stop("Argument 'n_blocks' must be defined.")
  }
  if (missing(initial_distance)) {
    stop("Argument 'initial_distance' must be defined.")
  }
  if (missing(increase)) {
    stop("Argument 'increase' must be defined.")
  }
  if (max_n_samplings > replicates) {
    stop("Argument 'replicates' must be larger than 'max_n_samplings'.")
  }
  if (is.null(master$master_matrix$Block)) {
    stop("Blocks are not defined in 'master', see function 'make_blocks'.")
  }
  if (!select_point[1] %in% c("random", "E_centroid", "G_centroid")) {
    stop("Argument 'select_point' is not valid.")
  }
  if (!cluster_method[1] %in% c("hierarchical", "k-means")) {
    stop("Argument 'cluster_method' is not valid.")
  }

  # running in loop for multiple answers if needed
  all_sites <- lapply(1:max_n_samplings, function(x) {
    ss <- set_seed + x - 1

    rule <- block_sample(master, variable_1, variable_2, n_blocks,
                         selection_type = "uniform", initial_distance, increase,
                         replicates, set_seed = ss)$master_matrix$Selected_blocks
    rule <- rule == 1

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
    rbind(unselp, ueselp, meselp)
  })

  names(all_sites) <- paste0("selection_", 1:length(all_sites))

  master$selected_sites_EG <- all_sites

  return(structure(master, class = "master_selection"))
}
