#' biosurvey: Tools for Biological Survey Planning
#'
#' biosurvey is a collection of tools that allow users to plan systems of sampling
#' sites increasing the efficiency of biodiversity monitoring by considering the
#' relationship between environmental and geographic conditions in a region.
#' Three main modules are included: 1) Data preparation; 2) Selection of sets
#' of sites for biodiversity sampling; and, 3) Tools for testing efficiency of
#' distinct sets of sampling sites. Data are prepared ways that avoid the need
#' for more data in posterior analyses, and allow concentrating in critical
#' methodological decisions to select sampling sites. Various algorithms for
#' selecting sampling sites are available, and options for considering
#' pre-selected sites (known to be important for biodiversity monitoring) are
#' included. Visualization is a critical component in this set of tools and most
#' of the results obtained can be plotted to help to understand their implications.
#' The options for selecting sampling sites include here differ from other
#' implementations in that they consider the environmental and geographic
#' structure of a region to suggest sampling sites that could increase the
#' efficiency of efforts dedicated to monitoring biodiversity.
#'
#' @section Main functions in biosurvey:
#' \code{\link{base_PAM}}, \code{\link{block_sample}},
#' \code{\link{compare_SAC}}, \code{\link{EG_selection}},
#' \code{\link{explore_data_EG}}, \code{\link{make_blocks}},
#' \code{\link{master_matrix}}, \code{\link{PAM_indices}},
#' \code{\link{plot_blocks_EG}}, \code{\link{plot_SAC}},
#' \code{\link{plot_sites_EG}}, \code{\link{random_selection}},
#' \code{\link{selected_sites_SAC}}, \code{\link{subset_PAM}},
#' \code{\link{uniformE_selection}}, \code{\link{uniformG_selection}}
#'
#'
#' Other functions (important helpers)
#'
#' \code{\link{assign_blocks}}, \code{\link{closest_to_centroid}},
#' \code{\link{distance_filter}}, \code{\link{files_2data}},
#' \code{\link{find_clusters}}, \code{\link{find_modes}},
#' \code{\link{grid_from_region}}, \code{\link{match_rformat}},
#' \code{\link{PAM_from_table}}, \code{\link{point_sample}},
#' \code{\link{point_sample_cluster}}, \code{\link{point_thinning}},
#' \code{\link{refill_PAM_indices}}, \code{\link{rlist_2data}},
#' \code{\link{selected_sites_PAM}}, \code{\link{spdf_2data}},
#' \code{\link{stack_2data}}, \code{\link{unimodal_test}},
#' \code{\link{wgs84_2aed_laea}}
#'
#' @docType package
#' @name biosurvey
NULL
