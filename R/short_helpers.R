#' Test wheather a number is pair
#'
#' @param x (numeric) value to be tested.
#' @return Logical value
#' @export
#' @examples
#' is_pair(4)
#' is_pair(5)

is_pair <- function(x) {
  x / 2 == as.integer(x / 2)
}
