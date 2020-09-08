context("Preparation of PAM subset")

test_that("Correct creation of PAM_subset", {
  sub_pam_G <- subset_PAM(b_pam, m_selection, selection_type = "G")
  sub_pam_all <- subset_PAM(b_pam, m_selection, selection_type = "all")

  cnam <- names(sub_pam_G)

  enams <- c("PAM", "PAM_indices", "PAM_selected_sites_random",
             "PAM_selected_sites_G", "PAM_selected_sites_E",
             "PAM_selected_sites_EG")

  testthat::expect_s3_class(sub_pam_G, "PAM_subset")
  testthat::expect_s4_class(sub_pam_G$PAM, "SpatialPolygonsDataFrame")
  testthat::expect_length(sub_pam_G, 6)
  testthat::expect_equal(cnam, enams)
  testthat::expect_null(sub_pam_G$PAM_selected_sites_random)
  testthat::expect_null(sub_pam_G$PAM_selected_sites_E)
  testthat::expect_null(sub_pam_G$PAM_selected_sites_EG)
  testthat::expect_null(sub_pam_all$PAM_selected_sites_EG)
})


test_that("Errors PAM_subset", {
  testthat::expect_error(subset_PAM(b_pam))
  testthat::expect_error(subset_PAM(master_selection = m_selection))
  testthat::expect_error(subset_PAM(b_pam, m_selection, selection_type = "i"))
})
