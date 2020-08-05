#' Selection of survey sites maximizing uniformity in environmental space
#' considering geographic structure
#'
#' @description Selection of sites to be sampled in a survey, with the goal of
#' maximizing uniformity of points in environment, but considering geographic
#' patterns of data. Sets of points that are environmentally similar and have a
#' disjoint pattern in geography, are selected twice (two survey sites are placed
#' so they consider the biggest geographic).
#'
#' @param master a master_matrix object derived from the function
#' \code{\link{master_matrix}} or a master_selection object derived from functions
#' \code{\link{random_selection}}, \code{\link{uniformG_selection}},
#' or \code{\link{uniformE_selection}}.
#' @param variable_1 (character or numeric) name or position of the first
#' variable (X-axis).
#' @param variable_2 (character or numeric) name or position of the second
#' variable (Y-axis).
#' @param n_blocks (numeric) number of blocks to be selected to be used as the
#' base for further explorations. Default = NULL.
#' @param initial_distance (numeric) euclidean distance to be used for a first
#' process of thinning and detection of remaining blocks. See details in
#' \code{\link{point_thinning}}.
#' @param increase (numeric) value to be added to \code{initial_distance} until
#' reaching the number in \code{n_blocks}.
#' @param replicates (numeric) number of thinning replicates performed to select
#' blocks uniformly. Default = 10.
#' @param use_preselected_sites (logical) whether to use sites that have been
#' defined as part of the selected sites previous any selection. Object in
#' \code{master} must contain the site(s) preselected in and element of name
#' "preselected_sites" for this argument to be effective. Default = TRUE.
#' See details for more information on the approach used.
#' @param n_optimization (numeric) number of times the algorithm of optimization
#' for site selection when using \code{use_preselected_sites} will run. Default
#' = 1000.
#' @param max_n_samplings (numeric) maximum number of samples to be chosen after
#' performing all thinning \code{replicates}. Default = 1.
#' @param select_point (character) how or which point will be selected. Three
#' options are available: "random", "E_centroid", "G_centroid". E_ or G_ centroid
#' indicate that the point(s) closets to the respective centroid will be selected.
#' Default = "E_centroid".
#' @param cluster_method (character) name of the method to be used for detecting
#' clusters of block points in geographic space. Options are "hierarchical" and
#' "k-means"; default = "hierarchical". See details in \code{\link{find_clusters}}.
#' @param sample_for_distance (numeric) sample to be considered when measuring
#' the geographic distances among points in the blocks of environmental points.
#' The distances measured are then used to test whether points are distributed
#' uniformly or not in the geography. Default = 250.
#' @param median_distance_filter (character) optional argument to define a median
#' distance-based filter based on which sets of sampling sites will be selected.
#' The default, NULL, does not apply such a filter. Options are: "max" and "min".
#' See details.
#' @param set_seed (numeric) integer value to specify a initial seed. Default = 1.
#'
#' @return
#' A master_selection object (S3) with a special element called
#' selected_sites_EG containing one or more sets of selected sites depending on
#' \code{max_n_samplings} and \code{median_distance_filter}.
#'
#' @details
#' Two important steps are needed before using this function: 1) exploring data in
#' environmental and geographic spaces, and 2) performing a regionalization of the
#' environmental space. Exploring the data can be done using the function
#' \code{\link{explore_data_EG}}. This step is optional but strongly recommended,
#' as important decisions that need to be taken depend on the configuration
#' of the data in the two spaces. A regionalization of the environmental space
#' of the region of interest helps in defining important parts of your region
#' that should be considered to select sites. This can be done using the function
#' \code{\link{make_blocks}}. Later the regions created in environmental space
#' will be used for selecting one or more sampling sites per block depending on
#' the geographic pattern of such environmental combinations.
#'
#' The process of survey-site selection with this function is the most complex
#' among all functions in this package. The complexity derives from the aim of the
#' function, which is to select sites that sample appropriately environmental
#' combinations in the region of interest (environmental space), but
#' considering the geographic patterns of such environmental regions (geographic
#' space).
#'
#' In this approach, the first step is to select candidate blocks (from the
#' ones obtained with \code{\link{make_blocks}}) that are uniformly distributed
#' in environmental space. The geographic configuration of points in such
#' blocks is explored to detect whether they are clustered (i.e., similar
#' environmental conditions are present in distant places in the region of
#' interest). For blocks with points that present one cluster in geography,
#' only one survey site is selected, and for those with multiple clusters in
#' geographic space, two survey sites are selected considering the two largest
#' clusters.
#'
#' As multiple sets could result from selection, the argument of the function
#' \code{median_distance_filter} could be used to select the set of sites with
#' the maximum ("max") or minimum ("min") median distance among selected sites.
#' Option "max" will increase the geographic distance among sampling sites, which
#' could be desirable if the goal is to cover the region of interest more broadly.
#' The other option "min", could be used in cases when the goal is to reduce
#' resources and time needed to sample such sites.
#'
#' @seealso
#' \code{\link{uniformG_selection}}, \code{\link{uniformE_selection}},
#' \code{\link{random_selection}}, \code{\link{make_blocks}},
#' \code{\link{plot_sites_EG}}
#'
#' @usage
#' EG_selection(master, variable_1, variable_2, n_blocks, initial_distance,
#'              increase, replicates = 10, use_preselected_sites = TRUE,
#'              n_optimization = 1000, max_n_samplings = 1,
#'              select_point = "E_centroid", cluster_method = "hierarchical",
#'              median_distance_filter = NULL, sample_for_distance = 250,
#'              set_seed = 1)
#'
#' @export
#'
#' @examples
#' \donttest{
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
#' }


EG_selection <- function(master, variable_1, variable_2, n_blocks,
                         initial_distance, increase, replicates = 10,
                         use_preselected_sites = TRUE, n_optimization = 1000,
                         max_n_samplings = 1, select_point = "E_centroid",
                         cluster_method = "hierarchical",
                         median_distance_filter = NULL,
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
  if (use_preselected_sites == TRUE & is.null(master$preselected_sites)) {
    message("Element 'preselected_sites' in 'master' is NULL, setting\n'use_preselected_sites' = FALSE.")
    use_preselected_sites <- FALSE
  }
  if (use_preselected_sites == FALSE) {
    if (missing(initial_distance)) {
    stop("Argument 'initial_distance' must be defined.")
    }
    if (missing(increase)) {
      stop("Argument 'increase' must be defined.")
    }
    if (max_n_samplings > replicates) {
      stop("Argument 'replicates' must be larger than 'max_n_samplings'.")
    }
  }
  if (use_preselected_sites == TRUE & is.null(master$preselected_sites$Block)) {
    stop("Blocks are not defined in 'preselected_sites', see function 'make_blocks'.")
  }
  if (is.null(master$master_matrix$Block)) {
    stop("Blocks are not defined in 'master_matrix', see function 'make_blocks'.")
  }
  if (!select_point[1] %in% c("random", "E_centroid", "G_centroid")) {
    stop("Argument 'select_point' is not valid, see function's help.")
  }
  if (!cluster_method[1] %in% c("hierarchical", "k-means")) {
    stop("Argument 'cluster_method' is not valid, see function's help.")
  }
  if (!is.null(median_distance_filter)) {
    if (!median_distance_filter %in% c("max", "min")) {
      stop("Argument 'median_distance_filter' is not valid, see function's help.")
    }
  }

  # running
  if (use_preselected_sites == TRUE){
    # using preselected sites, optimization based on random selection
    pre <- master$preselected_sites

    ## blocks preselected sites belong to
    preblock <- unique(pre$Block)
    blocks <- unique(master$master_matrix$Block)
    lim <- (length(preblock) + 1):n_blocks
    nlim <- length((length(preblock) + 1):n_blocks)

    ## finding most uniform randomly selected blocks
    sblocks <- lapply(1:n_optimization, function(x) {
      ss <- set_seed + x - 1
      selb <- sample(blocks, n_blocks)
      ub <- unique(c(preblock, selb))[lim]
    })

    ## keeping only unique sets
    cd <- sapply(sblocks, function(x) {paste0(sort(x), collapse = "_")})

    sblocks <- sblocks[which(!duplicated(cd))]

    ## obtaining block centroids (closest to centroid) info
    sblocks <- lapply(sblocks, function(x) {
      closest_to_centroid(master$master_matrix[master$master_matrix$Block ==
                                                 x, ], variable_1, variable_2,
                          space = "E", n = 1, id_column = "Block")
    })

    ## keeping only sblocks with correct number of blocks
    sblocks[sapply(sblocks, is.null)] <- NULL
    cd <- unlist(sapply(sblocks, nrow))
    sblocks <- sblocks[which(cd == nlim)]

    # Post filtering to get higher uniformity of blocks in E
    sblocks <- distance_filter(sblocks, "max")[[1]]
    sblocks <- unique(sblocks$Block)
    rules <- 1
  } else {
    # no user sites
    ## creating rule for block selection
    rules <- lapply(1:replicates, function(x) {
      ss <- set_seed + x - 1
      rule <- suppressMessages(block_sample(master, variable_1, variable_2, n_blocks,
                                            selection_type = "uniform",
                                            initial_distance, increase, replicates = 1,
                                            set_seed = ss)$master_matrix$Selected_blocks)
      which(rule == 1)
    })

    ## keeping only unique sets
    cd <- sapply(rules, function(x) {paste0(sort(x), collapse = "_")})

    rules <- rules[which(!duplicated(cd))]
  }

  # defining columns with coordinates
  g_cols <- c("Longitude", "Latitude")

  # finding sites
  all_sites <- lapply(rules, function(x) {
    ## subsetting based on rules
    if (use_preselected_sites == TRUE) {
      blocksp <- sblocks
    } else {
      blocksp <- unique(master$master_matrix[x, "Block"])
    }

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
    if (use_preselected_sites == TRUE) {
      rbind(pre[, -1], unselp, ueselp, meselp)
    } else {
      rbind(unselp, ueselp, meselp)
    }
  })

  # Getting needed samples form the most numerous ones
  lns <- sapply(all_sites, nrow)
  all_sites <- all_sites[which(lns == max(lns))]

  cd <- sapply(all_sites, function(x) {
    y <- x[order(x[, variable_1], x[, variable_2]), ]
    paste0(paste0(y[, variable_1], y[, variable_2]), collapse = "_")
  })

  ## excluding duplicates
  all_sites <- all_sites[which(!duplicated(cd))]
  nsel <- ifelse(length(all_sites) < max_n_samplings, length(all_sites),
                 max_n_samplings)

  # results to be returned
  if (nsel == 1) {
    all_sites <- list(all_sites[[1]])
  } else {
    all_sites <- all_sites[1:nsel]
  }

  # Select final set based on median geographic distance
  if (length(all_sites) > 1 & !is.null(median_distance_filter)) {
    all_sites <- distance_filter(all_sites, median_distance_filter)
  }

  # Preparing and returning results
  names(all_sites) <- paste0("selection_", 1:length(all_sites))

  master$selected_sites_EG <- all_sites

  return(structure(master, class = "master_selection"))
}
