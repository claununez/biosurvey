#' Simple color palettes
#'
#' @param n number of colors to be in the palette.
#'
#' @return
#' A character vector of hex color codes.
#'
#' @aliases purplow darkros daright greeple bluered
#'
#' @usage
#' # sequential palettes
#' purplow(n)
#'
#' @rdname purplow
#'
#' @examples
#' purplow(3)
#' darkros(3)
#' daright(3)
#' greeple(3)
#' bluered(3)

purplow <- colorRampPalette(rev(c("#edf8b1", "#41b6c4", "#081d58")))


#' @rdname purplow
#' @usage
#' darkros(n)
darkros <- colorRampPalette(rev(c("#fff7f3", "#f768a1", "#49006a")))


#' @rdname purplow
#' @usage
#' daright(n)
daright <- colorRampPalette(rev(c("#e0ecf4", "#8c96c6", "#4d004b")))


#' @rdname purplow
#' @usage
#' # diverging palettes
#' greeple(n)
greeple <- colorRampPalette(rev(c("#762a83", "#f7f7f7", "#1b7837")))


#' @rdname purplow
#' @usage
#' bluered(n)
bluered <- colorRampPalette(rev(c("#d73027", "#ffffbf", "#4575b4")))
