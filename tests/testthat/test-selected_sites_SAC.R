context("Preparation of SAC for sites")

test_that("Correct creation of SACs for sites", {
  sub_pam_all <- subset_PAM(b_pam, m_selection, selection_type = "all")
  SACs <- selected_sites_SAC(PAM_subset = sub_pam_all, selection_type = "all")

  cnam <- names(SACs)

  enams <- c("SAC_selected_sites_random", "SAC_selected_sites_E",
             "SAC_selected_sites_G")

  testthat::expect_length(SACs, 3)
  testthat::expect_equal(cnam, enams)
  testthat::expect_s3_class(SACs$SAC_selected_sites_G$selection_1, "specaccum")
  testthat::expect_s3_class(SACs$SAC_selected_sites_E$selection_1, "specaccum")
  testthat::expect_s3_class(SACs$SAC_selected_sites_random$selection_1, "specaccum")
})


test_that("Errors SACs for sites", {
  testthat::expect_error(selected_sites_SAC())
  testthat::expect_error(selected_sites_SAC(PAM_subset = sub_pam_all,
                                            selection_type = "i"))
})
