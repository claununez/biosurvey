#' Selection of survey sites maximizing uniformity in environmental space
#' considering geographic structure
#'
#' @description Selection of sites to be sampled in a survey, with the goal of
#' maximizing uniformity of points in the environment, but considering
#' geographic patterns of data. Sets of points that are environmentally similar
#' and have a disjoint pattern in geography, are selected twice (two survey
#' sites are placed so they consider the biggest geographic clusters).
#'
#' @param master master_matrix object derived from the function
#' \code{\link{prepare_master_matrix}} or master_selection object derived from
#' functions \code{\link{random_selection}}, \code{\link{uniformG_selection}},
#' or \code{\link{uniformE_selection}}.
#' @param n_blocks (numeric) number of blocks to be selected to be used as the
#' base for further explorations. If preselected sites are used, this number must
#' be larger than the number of unique blocks already represented by such sites.
#' @param guess_distances (logical) whether or not to use internal algorithm
#' to automatically select \code{initial_distance} and \code{increase}. Default
#' = TRUE. If FALSE, \code{initial_distance} and \code{increase} must be
#' defined.
#' @param initial_distance (numeric) Euclidean distance to be used for a first
#' process of thinning and detection of remaining blocks. See details in
#' \code{\link{point_thinning}}. Default = NULL.
#' @param increase (numeric) initial value to be added to or subtracted from
#' \code{initial_distance} until reaching the number of \code{expected_points}.
#' Default = NULL.
#' @param max_n_samplings (numeric) maximum number of samples to be chosen
#' after performing all thinning \code{replicates}. Default = 1.
#' @param replicates (numeric) number of thinning replicates performed to
#' select blocks uniformly. Default = 10.
#' @param use_preselected_sites (logical) whether to use sites that have been
#' defined as part of the selected sites previous any selection. Object in
#' \code{master} must contain the site(s) preselected in and element of name
#' "preselected_sites" for this argument to be effective. Default = TRUE.
#' See details for more information on the approach used.
#' @param select_point (character) how or which point will be selected for each
#' block or cluster. Three options are available: "random", "E_centroid", and
#' "G_centroid". E_ or G_ centroid indicate that the point(s) closets to the
#' respective centroid will be selected. Default = "E_centroid".
#' @param cluster_method (character) name of the method to be used for detecting
#' geographic clusters of points inside each block. Options are "hierarchical"
#' and "k-means"; default = "hierarchical". See details in
#' \code{\link{find_clusters}}.
#' @param sample_for_distance (numeric) sample to be considered when measuring
#' the geographic distances among points in blocks created in environmental
#' space. The distances measured are then used to test whether points are
#' distributed uniformly or not in the geography. Default = 250.
#' @param median_distance_filter (character) optional argument to define a
#' median distance-based filter based on which sets of sampling sites will be
#' selected. The default, NULL, does not apply such a filter. Options are:
#' "max" and "min". See details.
#' @param set_seed (numeric) integer value to specify a initial seed.
#' Default = 1.
#' @param verbose (logical) whether or not to print messages about the process.
#' Default = TRUE.
#' @param force (logical) whether to replace existing set of sites selected
#' with this method in \code{master}.
#'
#' @return
#' A \code{\link{master_selection}} object (S3) with a special element called
#' selected_sites_EG containing one or more sets of selected sites depending on
#' \code{max_n_samplings} and \code{median_distance_filter}.
#'
#' @details
#' Two important steps are needed before using this function: 1) exploring data
#' in environmental and geographic spaces, and 2) performing a regionalization
#' of the environmental space. Exploring the data can be done using the function
#' \code{\link{explore_data_EG}}. This step is optional but strongly
#' recommended, as important decisions that need to be taken depend on the
#' of the data in the two spaces. A regionalization of the environmental space
#' configuration of the region of interest helps in defining important parts of
#' your region that should be considered to select sites. This can be done
#' using the function \code{\link{make_blocks}}. Later, the regions created in
#' environmental space will be used for selecting one or more sampling sites per
#' block depending on the geographic pattern of such environmental combinations.
#'
#' The process of survey-site selection with this function is the most complex
#' among all functions in this package. The complexity derives from the aim of
#' the function, which is to select sites that sample appropriately
#' environmental combinations in the region of interest (environmental space),
#' but considering the geographic patterns of such environmental regions
#' (geographic space).
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
#' If \code{use_preselected_sites} is TRUE and such sites are included as an
#' element in the object in \code{master}, the approach for selecting sites in
#' environmental space considering geographic patterns is a little  different.
#' User-preselected sites will always be part of the sites selected. Other points
#' are selected based on an algorithm that searches for sites that are uniformly
#' distributed in environmental space but at a distance from preselected sites
#' that helps in maintaining uniformity among environmental blocks selected.
#' Note that preselected sites will not be processed, therefore, uniformity of
#' blocks representing such points cannot be warrantied.
#'
#' As multiple sets could result from selection, the argument of the function
#' \code{median_distance_filter} could be used to select the set of sites with
#' the maximum ("max") or minimum ("min") median distance among selected sites.
#' Option "max" will increase the geographic distance among sampling sites,
#' which could be desirable if the goal is to cover the region of interest more
#' broadly. The other option, "min", could be used in cases when the goal is to
#' reduce resources and time needed to sample such sites.
#'
#' @seealso
#' \code{\link{uniformG_selection}}, \code{\link{uniformE_selection}},
#' \code{\link{random_selection}}, \code{\link{make_blocks}},
#' \code{\link{plot_sites_EG}}
#'
#' @usage
#' EG_selection(master, n_blocks, guess_distances = TRUE, initial_distance = NULL,
#'              increase = NULL, max_n_samplings = 1, replicates = 10,
#'              use_preselected_sites = TRUE, select_point = "E_centroid",
#'              cluster_method = "hierarchical", median_distance_filter = NULL,
#'              sample_for_distance = 250, set_seed = 1,
#'              verbose = TRUE, force = FALSE)
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
#' colnames(m_blocks$data_matrix)
#'
#' # Selecting sites uniformly in E and G spaces
#' EG_sel <- EG_selection(master = m_blocks, n_blocks = 10,
#'                        initial_distance = 1.5, increase = 0.1,
#'                        replicates = 1, max_n_samplings = 1,
#'                        select_point = "E_centroid",
#'                        cluster_method = "hierarchical",
#'                        sample_for_distance = 100)
#'
#' head(EG_sel$selected_sites_EG[[1]])
#' dim(EG_sel$selected_sites_EG[[1]])
#' }


EG_selection <- function(master, n_blocks, guess_distances = TRUE,
                         initial_distance = NULL, increase = NULL,
                         max_n_samplings = 1, replicates = 10,
                         use_preselected_sites = TRUE,
                         select_point = "E_centroid",
                         cluster_method = "hierarchical",
                         median_distance_filter = NULL,
                         sample_for_distance = 250,
                         set_seed = 1, verbose = TRUE, force = FALSE) {

  # Initial tests
  if (missing(master)) {
    stop("Argument 'master' must be defined")
  }
  clsm <- class(master)[1]
  if (clsm %in% c("master_matrix", "master_selection")) {
    if (clsm == "master_selection") {
      if (!is.null(master$selected_sites_EG) & force == FALSE) {
        stop("'master' already contains a selection of this type, use 'force' = TRUE to replace it")
      }
    }
  } else {
    stop("Argument 'master' must be of class 'master_matrix' or 'master_selection'")
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
  if (missing(n_blocks)) {
    stop("Argument 'n_blocks' must be defined.")
  }
  if (use_preselected_sites == TRUE & is.null(master$preselected_sites)) {
    if (verbose == TRUE) {
      message("Element 'preselected_sites' in 'master' is NULL, setting\n'use_preselected_sites' = FALSE")
    }
    use_preselected_sites <- FALSE
  }
  if (guess_distances == FALSE) {
    if (missing(initial_distance)) {
      stop("Argument 'initial_distance' is not defined.")
    }
    if (missing(increase)) {
      stop("Argument 'increase' is not defined.")
    }
  }
  if (max_n_samplings > replicates) {
    stop("Argument 'replicates' must be larger than 'max_n_samplings'.")
  }
  if (use_preselected_sites == TRUE & is.null(master$preselected_sites$Block)) {
    stop("Blocks are not defined in 'preselected_sites', see function 'make_blocks'.")
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

  # Running
  if (verbose == TRUE) {
    message("Preparing data for analysis")
  }

  if (use_preselected_sites == TRUE) {
    # Using preselected sites to create mask and define distance
    tst <- preselected_dist_mask(master, expected_points = n_blocks,
                                 space = "E", variable_1 = variable_1,
                                 variable_2 = variable_2, use_blocks = TRUE,
                                 verbose = verbose)

    # Excluding points in close blocks from analysis
    cents <- closest_to_centroid(master$data_matrix, variable_1, variable_2,
                                 space = "E", n = 1, id_column = "Block")
    data <- sp::SpatialPointsDataFrame(cents[, c(variable_1, variable_2)], cents,
                                       proj4string = sp::CRS("+init=epsg:4326"))
    bqs <- data[is.na(sp::over(data, tst$mask)$ID), ]@data$Block

    ## Preparing temporal master matrix to be used in new selection
    data <- master$data_matrix
    data <- data[data$Block %in% bqs, ]
    tmm <- new_master_matrix(data_matrix = data, region = master$region,
                             raster_base = master$raster_base)

    # npre
    pre <- master$preselected_sites
    n_blocks <- n_blocks - length(unique(pre$Block))
  }

  #  Creating rule for block selection
  if (verbose == TRUE) {
    message("Selecting relevant environmental blocks, please wait...")
  }
  rules <- lapply(1:replicates, function(x) {
    ss <- set_seed + x - 1
    if (use_preselected_sites == TRUE) {
      rule <- block_sample(tmm, n_blocks,
                           selection_type = "uniform", replicates = 1,
                           set_seed = ss)$data_matrix
    } else {
      rule <- block_sample(master, n_blocks,
                           selection_type = "uniform", replicates = 1,
                           set_seed = ss)$data_matrix
    }
    unique(rule[rule$Selected_blocks == 1, "Block"])
  })

  ## Keeping only unique sets
  cd <- sapply(rules, function(x) {paste0(sort(x), collapse = "_")})

  rules <- rules[which(!duplicated(cd))]

  # Defining columns with coordinates
  g_cols <- c("Longitude", "Latitude")

  # Finding sites
  if (verbose == TRUE) {
    message("Running algorithm for selecting sites, please wait...")
  }
  all_sites <- lapply(1:length(rules), function(x) {
    ## Measuring distances
    distsp <- lapply(rules[[x]], function(y) {
      block_data <- master$data_matrix[master$data_matrix[, "Block"] == y,
                                         g_cols]
      tpoints <- nrow(block_data)

      if (tpoints > sample_for_distance) {
        block_data <- block_data[sample(tpoints, sample_for_distance), ]
      }

      dsnna <- na.omit(c(raster::pointDistance(block_data, lonlat = TRUE)))
      dsnna[dsnna != 0]
    })
    names(distsp) <- rules[[x]]

    ## Unimodal tests and splitting data according to results
    dpp <- unimodal_test(distsp)

    unimp <- dpp[which(dpp$p_alue > 0.05), ]
    mmodp <- dpp[which(dpp$p_alue <= 0.05), ]
    nmodp <- dpp[which(is.na(dpp$p_alue)), ]

    ## Analysis with no mode (very few points)
    if (nrow(nmodp) > 0) {
      unselp <- point_sample(master$data_matrix[master$data_matrix[, "Block"]
                                                %in% nmodp[, "Block"], ],
                             variable_1, variable_2, n = 1,
                             select_point = "random", id_column = "Block")
    } else {
      unselp <- matrix(nrow = 0, ncol = ncol(master$data_matrix))
      colnames(unselp) <- colnames(master$data_matrix)
      unselp <- as.data.frame(unselp)
    }

    ## Analysis with unimodal
    if (nrow(unimp) > 0) {
      ueselp <- point_sample(master$data_matrix[master$data_matrix[, "Block"]
                                                %in%  unimp[, "Block"], ],
                             variable_1, variable_2, n = 1,
                             select_point = select_point, id_column = "Block")
    } else {
      ueselp <- matrix(nrow = 1, ncol = ncol(master$data_matrix))
      colnames(ueselp) <- colnames(master$data_matrix)
      ueselp <- as.data.frame(na.omit(ueselp))
    }

    ## Analysis with multimodal
    if (nrow(mmodp) > 0) {
      distsp <- distsp[as.character(mmodp$Block)]
      meselp <- point_sample_cluster(master$data_matrix[master$data_matrix[, "Block"]
                                                       %in% mmodp[, "Block"], ],
                                     variable_1, variable_2,
                                     distance_list = distsp,
                                     n = 1, cluster_method = cluster_method,
                                     select_point = select_point,
                                     id_column = "Block")
    } else {
      meselp <- matrix(nrow = 0, ncol = ncol(master$data_matrix))
      colnames(meselp) <- colnames(master$data_matrix)
      meselp <- as.data.frame(meselp)
    }

    if (verbose == TRUE) {
      message("    Process ", x, " of ", length(rules))
    }

    ## Combining all results
    if (use_preselected_sites == TRUE) {
      rbind(pre[, -1], unselp[, colnames(unselp) != "Selected_blocks"],
            ueselp[, colnames(ueselp) != "Selected_blocks"],
            meselp[, colnames(meselp) != "Selected_blocks"])
    } else {
      rbind(unselp[, colnames(unselp) != "Selected_blocks"],
            ueselp[, colnames(ueselp) != "Selected_blocks"],
            meselp[, colnames(meselp) != "Selected_blocks"])
    }
  })

  # Keeping the most numerous ones
  lns <- sapply(all_sites, nrow)
  all_sites <- all_sites[which(lns == max(lns))]

  ## Excluding duplicates
  cd <- sapply(all_sites, function(x) {
    y <- x[order(x[, variable_1], x[, variable_2]), ]
    paste0(paste0(y[, variable_1], y[, variable_2]), collapse = "_")
  })

  all_sites <- all_sites[which(!duplicated(cd))]

  # Results to be returned
  nsel <- ifelse(length(all_sites) < max_n_samplings, length(all_sites),
                 max_n_samplings)

  if (nsel == 1) {
    all_sites <- list(all_sites[[1]])
  } else {
    all_sites <- all_sites[1:nsel]
  }

  # Select final set based on median geographic distance
  if (length(all_sites) > 1 & !is.null(median_distance_filter)) {
    all_sites <- distance_filter(all_sites, median_distance_filter)
  }

  if (verbose == TRUE) {
    message("Total number of sites selected: ", nrow(all_sites[[1]]))
  }


  # Preparing and returning results
  ## Naming and returning
  names(all_sites) <- paste0("selection_", 1:length(all_sites))

  ## Arguments as attributes
  other_args <- list(arguments = list(variable_1 = variable_1,
                                      variable_2 = variable_2,
                                      n_blocks = n_blocks,
                                      guess_distances = guess_distances,
                                      max_n_samplings = max_n_samplings,
                                      replicates = replicates,
                                      use_preselected_sites = use_preselected_sites,
                                      select_point = select_point,
                                      cluster_method = cluster_method,
                                      median_distance_filter = median_distance_filter,
                                      sample_for_distance = sample_for_distance,
                                      set_seed = set_seed))
  attributes(all_sites) <- c(attributes(all_sites), other_args)


  if (class(master)[1] == "master_matrix") {
    return(new_master_selection(master$data_matrix, master$preselected_sites,
                                master$region, master$mask, master$raster_base,
                                master$PCA_results, selected_sites_random = NULL,
                                selected_sites_G = NULL, selected_sites_E = NULL,
                                all_sites))

  } else {
    return(new_master_selection(master$data_matrix, master$preselected_sites,
                                master$region, master$mask, master$raster_base,
                                master$PCA_results, master$selected_sites_random,
                                master$selected_sites_G, master$selected_sites_E,
                                all_sites))
  }
}
