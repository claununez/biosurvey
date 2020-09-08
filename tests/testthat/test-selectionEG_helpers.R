#----
context("Sampling points")

test_that("Correct point_sample", {
  points_s <- point_sample(m_matrix$data_matrix, variable_1 = "Max_temperature",
                           variable_2 = "Min_temperature", n = 3,
                           select_point = "E_centroid", id_column = NULL)

  nr <- dim(points_s)
  cnam <- c(colnames(m_matrix$data_matrix), "id_column")

  testthat::expect_equal(nr, c(3, 11))
  testthat::expect_equal(cnam, colnames(points_s))
})


test_that("Errors point_sample", {
  testthat::expect_error(point_sample())
  testthat::expect_error(point_sample(m_matrix$data_matrix))
  testthat::expect_error(point_sample(m_matrix$data_matrix,
                                      variable_1 = "Max_temperature"))
  testthat::expect_error(point_sample(m_matrix$data_matrix,
                                      variable_2 = "Min_temperature"))
})
#----

#----
context("Sampling points G clustered")

test_that("Correct point_sample_cluster", {
  m_blocks <- make_blocks(m_matrix, variable_1 = "PC1", variable_2 = "PC2",
                          n_cols = 10, n_rows = 10, block_type = "equal_area")
  datam <- m_blocks$data_matrix
  datam <- datam[datam$Block %in% names(dist_list), ]
  point_clus <- point_sample_cluster(datam, variable_1 = "PC1", variable_2 = "PC2",
                                     distance_list = dist_list, n = 1,
                                     cluster_method = "hierarchical",
                                     select_point = "E_centroid", id_column = "Block")

  dtda <- dim(point_clus)
  cnam <- colnames(point_clus)

  testthat::expect_equal(dtda, c(9, 11))
  testthat::expect_equal(cnam, colnames(datam))
})


test_that("Errors point_sample_cluster", {
  m_blocks <- make_blocks(m_matrix, variable_1 = "PC1", variable_2 = "PC2",
                          n_cols = 10, n_rows = 10, block_type = "equal_area")
  datam <- m_blocks$data_matrix
  datam <- datam[datam$Block %in% names(dist_list), ]
  point_clus <- point_sample_cluster(datam, variable_1 = "PC1", variable_2 = "PC2",
                                     distance_list = dist_list, n = 1,
                                     cluster_method = "hierarchical",
                                     select_point = "E_centroid", id_column = "Block")
  testthat::expect_error(point_sample_cluster())
  testthat::expect_error(closest_to_centroid(datam))
  testthat::expect_error(closest_to_centroid(datam, variable_1 = "PC1"))
  testthat::expect_error(closest_to_centroid(datam, variable_1 = "PC1",
                                             variable_2 = "PC2"))
})
#----


#----
context("Find geofraphic clusters")

test_that("Correct find_clusters", {
  clusters <-  find_clusters(m_matrix$data_matrix, x_column = "PC1",
                             y_column = "PC2", space = "E",
                             cluster_method = "hierarchical", n_k_means = NULL,
                             split_distance = 4)

  dtda <- dim(clusters)
  cnam <- colnames(clusters)

  testthat::expect_equal(dtda, c(nrow(m_matrix$data_matrix), 11))
  testthat::expect_equal(cnam, c(colnames(m_matrix$data_matrix), "clusters"))
})


test_that("Errors find_clusters", {
  testthat::expect_error(find_clusters())
  testthat::expect_error(find_clusters(m_matrix$data_matrix))
  testthat::expect_error(find_clusters(m_matrix$data_matrix, x_column = "PC1"))
  testthat::expect_error(find_clusters(m_matrix$data_matrix, x_column = "PC1",
                                       y_column = "PC2"))
  testthat::expect_error(find_clusters(m_matrix$data_matrix, x_column = "PC1",
                                       y_column = "PC2", space = "E"))

})
#----
