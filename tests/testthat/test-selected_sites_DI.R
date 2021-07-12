context("Preparation of dissimilarity indices")


test_that("Correct creation of dissimilarity indices", {
  data("b_pam", package = "biosurvey")
  data("m_selection", package = "biosurvey")

  sub_pam_all <- subset_PAM(b_pam, m_selection, selection_type = "all")
  DI_sel <- selected_sites_DI(sub_pam_all)

  DI_class <- class(DI_sel)
  cla <- class(sub_pam_all)
  cnam <- names(DI_sel)
  anames <- c("DI_selected_sites_random", "cluster_random",
              "DI_selected_sites_G", "cluster_G", "DI_selected_sites_E",
              "cluster_E", "all_selections", "DI_selections",
              "cluster_selections")

  testthat::expect_equal(DI_class, "list")
  testthat::expect_equal(cla[2], "base_PAM")
  testthat::expect_length(DI_sel, 9)
  testthat::expect_equal(cnam, anames)
})


test_that("Errors and messages selected_sites_DI", {
  data("b_pam", package = "biosurvey")
  data("m_selection", package = "biosurvey")

  testthat::expect_error(selected_sites_DI())
  testthat::expect_error(selected_sites_DI(PAM_subset, selection_type = "i"))
})
