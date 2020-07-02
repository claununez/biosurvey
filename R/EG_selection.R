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
#' @param selection_option (character) option of site selection; options are:
#' "distance_rule" and "G_clusters". Default = "distance_rule". See details.
#' @param expected_points (numeric) if \code{selection_option} = "distance_rule",
#' number of survey points (sites) to be selected. Default = NULL.
#' @param n_blocks (numeric) if \code{selection_option} = "G_clusters",
#' number of blocks to be selected from all existent blocks in
#' \code{master$master_matrix}. Default = NULL.
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
#' See details in \code{\link{find_clusters}}.
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
#' @details
#' Two important steps are needed before using this function: exploring data in
#' environmental and geographic spaces, and performing a rationalization of the
#' environmental space. Exploring the data can be done using the function
#' \code{\link{explore_data_EG}}. This step is optional but strongly recommended,
#' as may important decisions that need to be taken depend on the configuration
#' of the data in the two spaces. A rationalization of the environmental space
#' of the region of interest helps in defining important parts of your region
#' that should be considered to select sites. This can be done using the function
#' \code{\link{make_blocks}}. Later the regions created in environmental space
#' will be used for selecting one or more sampling sites per block depending on
#' the geographic pattern of such environmental combinations.
#'
#' The process of survey-site selection with this function is the most complex
#' among all functions in this package. The complexity derives from the aim of the
#' function, which is to select sites that sample appropriately environmental
#' combinations in the region of interest (environmental space), but also
#' considering the geographic patterns of such environmental regions (geographic
#' space). Two options for selection are available:
#'
#' 1. "distance_rule".- In this option, multiple sets of sites are selected
#' aiming for uniform distributions in environmental space, then, geographic
#' distances are measured among points of candidate sets, and only the set or
#' sets with maximum median distances are kept.
#'
#' 2. "G_clusters"._ Here, the first step is to select candidate blocks (from the
#' ones obtained with \code{\link{make_blocks}}) that are uniformly distributed
#' in environmental space. The geographic configuration of points in such
#' blocks is explored to detect whether they are clustered (i.e., similar
#' environmental conditions are present in distant places in the region of
#' interest). For blocks with points that are not clustered in geographic space,
#' only one survey site is selected, and for those with clustered geographic
#' patterns, two survey sites are selected considering the largest clusters.
#'
#' @seealso
#' \code{\link{uniformG_selection}}, \code{\link{uniformE_selection}},
#' \code{\link{random_selection}}, \code{\link{make_blocks}},
#' \code{\link{plot_sites_EG}}
#'
#' @usage
#' EG_selection(master, variable_1, variable_2, selection_option = "distance_rule",
#'              n_blocks = NULL, initial_distance, increase, replicates = 10,
#'              max_n_samplings = 1, select_point = "E_centroid",
#'              cluster_method = "hierarchical", sample_for_distance = 250,
#'              set_seed = 1)
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


EG_selection <- function(master, variable_1, variable_2,
                         selection_option = "distance_rule",
                         expected_points = NULL, n_blocks = NULL,
                         initial_distance, increase, replicates = 10,
                         max_n_samplings = 1, select_point = "E_centroid",
                         cluster_method = "hierarchical",
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
  if (!selection_option %in% c("distance_rule", "G_clusters")) {
    stop("Argument 'selection_option' is not valid. See function's help.")
  }
  if (selection_option == "G_clusters") {
    if (is.null(n_blocks)) {
      stop("If 'selection_option' = 'G_clusters', argument 'n_blocks' must be defined.")
    }
  } else {
    if (is.null(expected_points)) {
      stop("If 'selection_option' = 'G_clusters', argument 'expected_points' must be defined.")
    }
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
    stop("Argument 'select_point' is not valid. See function's help.")
  }
  if (!cluster_method[1] %in% c("hierarchical", "k-means")) {
    stop("Argument 'cluster_method' is not valid.")
  }

  # running
  if (selection_option = "distance_rule") {
    all_sites <- uniformE_selection(master, variable_1, variable_2,
                                    selection_from = "all_points", expected_points,
                                    max_n_samplings, initial_distance, increase,
                                    replicates, set_seed)$selected_sites_E

    if (length(all_sites) > 1) {
      dists <- sapply(all_sites, function(x) {
        dis <- raster::pointDistance(x[, c("Longitude", "Latitude")], lonlat = TRUE)
        diag(dis) <- NA
        median(c(dis), na.rm = T)
      })
      all_sites <- all_sites[dists == max(dists)]
    }
  } else {
    ## creating rule for block selection
    rules <- lapply(1:replicates, function(x) {
      ss <- set_seed + x - 1
      rule <- suppressMessages(block_sample(master, variable_1, variable_2, n_blocks,
                                            selection_type = "uniform",
                                            initial_distance, increase, replicates = 1,
                                            set_seed = ss)$master_matrix$Selected_blocks)
      which(rule == 1)
    })

    cd <- sapply(rules, function(x) {paste0(sort(x), collapse = "_")})

    rules <- rules[which(!duplicated(cd))]

    g_cols <- c("Longitude", "Latitude") # defning columns with coordinates

    all_sites <- lapply(rules, function(x) {
      ## subsetting based on rule
      blocksp <- unique(master$master_matrix[x, "Block"])

      ## measuring distances
      distsp <- lapply(blocksp, function(y) {
        block_data <- master$master_matrix[master$master_matrix[, "Block"] == y,
                                           g_cols]
        tpoints <- nrow(block_data)

        if (tpoints > sample_for_distance) {
          block_data <- block_data[sample(tpoints, sample_for_distance), ]
        }

        dsnna <- na.omit(c(raster::pointDistance(block_data, lonlat = TRUE)))
        dsnna[dsnna != 0]
      })
      names(distsp) <- blocksp

      ## unimodal tests and splitting data according to results
      dpp <- unimodal_test(distsp)

      unimp <- dpp[which(dpp$p_alue > 0.05), ]
      mmodp <- dpp[which(dpp$p_alue <= 0.05), ]
      nmodp <- dpp[which(is.na(dpp$p_alue)), ]

      ## analysis with no mode (very few points)
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

      ## analysis with unimodal
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

      ## analysis with multimodal
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

      ## combining all results
      rbind(unselp, ueselp, meselp)
    })
  }

  names(all_sites) <- paste0("selection_", 1:length(all_sites))

  master$selected_sites_EG <- all_sites

  return(structure(master, class = "master_selection"))
}
