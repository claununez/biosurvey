#----
context("Explore data EG")

test_that("Errors explore_data_EG", {
  testthat::expect_error(explore_data_EG())
  testthat::expect_error(explore_data_EG(m_matrix))
  testthat::expect_error(explore_data_EG(m_matrix,
                                         variable_1 = "Mean_temperature"))
  testthat::expect_error(explore_data_EG(m_matrix,
                                         variable_2 = "Annual_precipitation"))
})
#----

#----
context("Plot blocks EG")

test_that("Errors plot_blocks_EG", {
  m_blocks <- make_blocks(m_matrix, variable_1 = "PC1",
                          variable_2 = "PC2", n_cols = 10, n_rows = 10,
                          block_type = "equal_area")

  testthat::expect_error(plot_blocks_EG())
  testthat::expect_error(plot_blocks_EG(master = m_matrix, variable_1 = "PC1",
                                        variable_2 = "PC2"))
  testthat::expect_error(plot_blocks_EG(master = m_blocks, variable_1 = "PC1"))
  testthat::expect_error(plot_blocks_EG(master = m_blocks, variable_2 = "PC2"))
})
#----


#----
context("Plot sites EG")

test_that("Errors plot_sites_EG", {
  m_blocks <- make_blocks(m_matrix, variable_1 = "PC1",
                          variable_2 = "PC2", n_cols = 10, n_rows = 10,
                          block_type = "equal_area")
  selectionE <- uniformE_selection(m_blocks, variable_1 = "PC1",
                                   variable_2 = "PC2",
                                   selection_from = "block_centroids",
                                   expected_points = 20, max_n_samplings = 1,
                                   replicates = 5)

  testthat::expect_error(plot_sites_EG())
  testthat::expect_error(plot_sites_EG(selectionE))
  testthat::expect_error(plot_sites_EG(selectionE, selection_type = "A"))
  testthat::expect_error(plot_sites_EG(selectionE, selection_type = "G"))
})
#----


#----
context("Plot SAC")

test_that("Errors plot_SAC", {
  sub_pam_all <- subset_PAM(b_pam, m_selection, selection_type = "all")
  SACs <- selected_sites_SAC(PAM_subset = sub_pam_all, selection_type = "all")

  testthat::expect_error(plot_SAC())
})
#----


#----
context("Plot compare SAC")

test_that("Errors compare_SAC", {
  sub_pam_all <- subset_PAM(b_pam, m_selection, selection_type = "all")
  SACs <- selected_sites_SAC(PAM_subset = sub_pam_all, selection_type = "all")

  testthat::expect_error(compare_SAC())
  testthat::expect_error(compare_SAC(SACs))
  testthat::expect_error(compare_SAC(SACs, element_1 = 1))
  testthat::expect_error(compare_SAC(SACs, element_2 = 2))
  testthat::expect_error(compare_SAC(SACs, element_1 = 1, element_2 = 4))
})
#----


#----
context("Plot PAM geography")

test_that("Errors plot_PAM_geo", {
  testthat::expect_error(plot_PAM_geo())
  testthat::expect_error(plot_PAM_geo(b_pam, index = "ZXY"))
})
#----


#----
context("Plot PAM CS diagram")

test_that("Errors plot_PAM_CS", {
  testthat::expect_error(plot_PAM_CS())
})
#----


#----
context ("Matrix-like plot of dissimilarities indices")

test_that("Errors plot_DI", {
  testthat::expect_error(plot_DI())
  testthat::expect_error(plot_DI(DI_selected_sites, selection_type = "i"))
})
#----


context ("Dendograme plot of dissimilarities indices")

test_that("Errors DI_dendrogram", {
  testthat::expect_error(DI_dendrogram())
  testthat::expect_error(DI_dendrogram(DI_selected_sites, selection_type = "i"))
})
#----
