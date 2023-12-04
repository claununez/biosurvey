#' Randomize matrix using the curve ball algorithm
#'
#' @description
#' Randomize a matrix using the curve ball algorithm. This keeps the sum of
#' rows and columns constant.
#'
#' @param matrix matrix to be randomized.
#'
#' @returns
#' A randomized matrix.
#'
#' @details
#' The curve ball algorithm was described by Strona et al.
#' (2014)[https://www.nature.com/articles/ncomms5114].

randomize_matrix_cb <- function(matrix) {
  RC <- dim(matrix)
  R <- RC[1]
  C <- RC[2]

  hp <- list()

  for (row in 1:dim(matrix)[1]) {
    hp[[row]] <- (which(matrix[row, ] == 1))
  }

  l_hp <- length(hp)

  for (rep in 1:(5 * l_hp)) {
    AB <- sample(1:l_hp, 2)
    a <- hp[[AB[1]]]
    b <- hp[[AB[2]]]
    ab <- intersect(a,b)
    l_ab <- length(ab)
    l_a <- length(a)
    l_b <- length(b)

    if ((l_ab %in% c(l_a, l_b)) == FALSE) {
      tot <- setdiff(c(a, b), ab)
      l_tot <- length(tot)
      tot <- sample(tot, l_tot, replace = FALSE, prob = NULL)
      L <- l_a - l_ab

      hp[[AB[1]]] <- c(ab, tot[1:L])
      hp[[AB[2]]] <- c(ab, tot[(L+1):l_tot])
    }
  }

  rm <- matrix(0, R, C)
  for (row in 1:R) {
    rm[row, hp[[row]]] <- 1
  }

  return(rm)
}
