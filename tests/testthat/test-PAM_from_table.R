context("PAM from table")

test_that("PAM created correctly", {
  pam <- PAM_from_table(data = sp_data, ID_column = "ID",
                        species_column = "Species")
  dpam <- dim(pam)
  cnam <- colnames(pam)
  uval <- unique(unlist(pam[, -1]))

  testthat::expect_equal(dpam, c(173, 26))
  testthat::expect_equal(cnam, c("ID", paste0("Species_", 1:25)))
  testthat::expect_equal(uval, c(0, 1))
})


test_that("Errors pft", {
  testthat::expect_error(PAM_from_table(data = sp_data, ID_column = "I",
                                        species_column = "Species"))
  testthat::expect_error(PAM_from_table(data = sp_data, ID_column = "ID",
                                        species_column = "Spes"))
  testthat::expect_error(PAM_from_table(ID_column = "ID",
                                        species_column = "Spes"))
  testthat::expect_error(PAM_from_table(data = sp_data,
                                        species_column = "Spes"))
  testthat::expect_error(PAM_from_table(data = sp_data, ID_column = "ID"))
})
