#' biosurvey: Tools for Biological Survey Planning
#'
#' @description
#' biosurvey is a collection of tools that allow users to plan sampling sites.
#' The methods presented increase the efficiency of biodiversity monitoring by
#' considering the relationship between environmental and geographic conditions
#' in a region.
#'
#' @details
#' Three main modules are included: 1) data preparation, 2) selection of sets
#' of sites for biodiversity sampling, and 3) tools for testing efficiency of
#' distinct sets of sampling sites. Data are prepared in ways that avoid the
#' need for more data in posterior analyses, and allow concentrating in critical
#' methodological decisions to select sampling sites. Various algorithms for
#' selecting sampling sites are available, and options for considering
#' pre-selected sites (known to be important for biodiversity monitoring) are
#' included. Visualization is a critical component in this set of tools and most
#' of the results obtained can be plotted to help to understand their
#' implications. The options for selecting sampling sites included here differ
#' from other implementations in that they consider the environmental and
#' geographic structure of a region to suggest sampling sites that could
#' increase the efficiency of efforts dedicated to monitoring biodiversity.
#'
#' @section Main functions in biosurvey:
#' \code{\link{block_sample}}, \code{\link{compare_SAC}},
#' \code{\link{EG_selection}}, \code{\link{DI_dendrogram}},
#' \code{\link{explore_data_EG}}, \code{\link{make_blocks}},
#' \code{\link{PAM_indices}}, \code{\link{plot_blocks_EG}},
#' \code{\link{plot_DI}}, \code{\link{plot_PAM_geo}},
#' \code{\link{plot_PAM_CS}}, \code{\link{plot_SAC}},
#' \code{\link{plot_sites_EG}}, \code{\link{prepare_base_PAM}},
#' \code{\link{prepare_master_matrix}}, \code{\link{prepare_PAM_CS}},
#' \code{\link{random_selection}}, \code{\link{read_master}},
#' \code{\link{read_PAM}}, \code{\link{save_master}}, \code{\link{save_PAM}},
#' \code{\link{selected_sites_DI}}, \code{\link{selected_sites_SAC}},
#' \code{\link{subset_PAM}}, \code{\link{uniformE_selection}},
#' \code{\link{uniformG_selection}}
#'
#' @section Other functions (important helpers):
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
#' @section Data included:
#' \code{\link{b_pam}}, \code{\link{distance_filter}}, \code{\link{dist_list}},
#' \code{\link{files_2data}}, \code{\link{m_matrix}},
#' \code{\link{m_matrix_pre}}, \code{\link{m_selection}}, \code{\link{mx}},
#' \code{\link{preselected}}, \code{\link{sp_data}}, \code{\link{species_data}},
#' \code{\link{sp_layers}}, \code{\link{sp_occurrences}},
#' \code{\link{variables}}
#'
#' @docType package
#' @name biosurvey
NULL
