#' Helper function to find raster extension
#'
#' @param format (character) any of the format names allowed for raster objects.
#' Options are: "GTiff", "ascii", "EHdr", "SAGA", "IDRISI", "CDF", "ENVI",
#' and "HFA".
#'
#' @return Raster extension according to format type.
#'
#' @export
#'
#' @examples
#' match_rformat("GTiff")

match_rformat <- function(format) {
  # Initial test
  if (missing(format)) {stop("Argument 'format' needs to be defined.")}
  # Defining format
  if (format == "GTiff") {format1 <- ".tif"}
  if (format == "EHdr") {format1 <- ".bil"}
  if (format == "ascii") {format1 <- ".asc"}
  if (format == "SAGA") {format1 <- ".sdat"}
  if (format == "IDRISI") {format1 <- ".rst"}
  if (format == "CDF") {format1 <- ".nc"}
  if (format == "ENVI") {format1 <- ".envi"}
  if (format == "HFA") {format1 <- ".img"}
  return(format1)
}





#' Project spatial points from geographic coordinates
#'
#' @param data matrix or data.frame that contains at least two columns, one
#' with longitude information and the other with latitude information.
#' @param longitude (character) the name of the column that contains the
#' longitude information.
#' @param latitude (character) the name of the column that contains the latitude
#' information.
#' @param which (character) type of projection. There are two options available:
#' "ED", for Azimuthal Equidistant and "EA", for Lambert Azimuthal Equal-Area.
#' Default = "ED".
#'
#' @return SpatVector projected to an option in \code{which}.
#'
#' @usage
#' wgs84_2aed_laea(data, longitude, latitude, which = "ED")
#'
#' @export
#' @importFrom terra crs vect project
#'
#' @examples
#' data("sp_occurrences", package = "biosurvey")
#'
#' sp_occ <- wgs84_2aed_laea(sp_occurrences, longitude = "longitude",
#'                           latitude = "latitude", which = "EA")

wgs84_2aed_laea <- function (data, longitude, latitude, which = "ED") {
  # Initial tests
  if (missing(data)) {
    stop("Argument 'data' is not defined.")
  }
  if (missing(longitude)) {
    stop("Argument 'longitude' is not defined.")
  }
  if (missing(latitude)) {
    stop("Argument 'latitude' is not defined.")
  }

  # Initial projection
  WGS84 <- terra::crs("+init=epsg:4326")
  dat_s <- terra::vect(data, geom = c(x = longitude, y = latitude),
                       crs = WGS84)

  # Projecting points
  cent <- apply(data[, c(longitude, latitude)], 2, mean)
  ini <- ifelse(which[1] == "ED", "+proj=aeqd", "+proj=laea")
  prj <- terra::crs(paste0(ini, " +lat_0=", cent[2], " +lon_0=", cent[1],
                        " +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

  dat_s <- terra::project(dat_s, prj)

  return(dat_s)
}




#' Helper to add a bar image legend to plots
#'
#' @param position (numeric or character) position of the bottom left corner of
#' the legend. If numeric, x and y coordinates. If character, options are:
#' "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright",
#' "right", or "center".
#' @param col color palette. A vector of contiguous colors. It can be generated
#' using functions like \code{\link{purplow}} (e.g., \code{darkros(255)}).
#' @param width_prop width of bar legend represented as a proportion of the
#' entire plotting width. Default = 0.03.
#' @param heigh_prop heigh of bar legend represented as a proportion of the
#' entire plotting heigh. Default = 0.18.
#' @param title legend title. Default = "Legend".
#' @param labels (numeric or character) labels for the legend. Default =
#' c("Low", "High").
#' @param digits (numeric) number of decimal places to round numeric labels.
#' Default = 0.
#' @param labels_offset offset of labels from bar. Default = 0.2.
#' @param horizontal (logical) should the legend be horizontal. Default = FALSE.
#' @param alpha (numeric) alpha level 0-1. Default = 1.
#' @param border color for the border of the legend bar. The default, NULL, does
#' not plot a border.
#' @param cex character expansion factor relative to current.
#' Default = NULL.
#' @param inset inset distances from plot margins relative to plot region.
#' Default = 0.05.
#' @param insetx inset from x margins. The default, NULL, uses \code{inset}.
#' @param insety inset from y margins. The default, NULL, uses \code{inset}.
#'
#' @return
#' A bar legend for a plot.
#'
#' @export
#'
#' @importFrom grDevices as.raster
#' @importFrom graphics rasterImage polygon
#'
#' @usage
#' legend_bar(position, col, width_prop = 0.03, heigh_prop = 0.18,
#'            title = "Legend", labels = c("Low", "High"), digits = 0,
#'            labels_offset = 0.2, horizontal = FALSE, alpha = 1, border = NULL,
#'            cex = NULL, inset = 0.05, insetx = NULL, insety = NULL)


legend_bar <- function(position, col, width_prop = 0.03, heigh_prop = 0.18,
                       title = "Legend", labels = c("Low", "High"), digits = 0,
                       labels_offset = 0.2, horizontal = FALSE, alpha = 1,
                       border = NULL, cex = NULL, inset = 0.05,
                       insetx = NULL, insety = NULL) {
  # Initial tests
  if (missing(position)) {
    stop("Argument 'position' is required to produce the legend")
  }
  if (missing(col)) {
    stop("Argument 'col' is required to produce the legend")
  }
  if (!class(position)[1] %in% c("numeric", "character")) {
    stop()
  } else {
    if (class(position)[1] == "character") {
      position <- match.arg(position, c("bottomright", "bottom", "bottomleft",
                                        "left", "topleft", "top", "topright",
                                        "right", "center"))
    }
  }

  # limits of plot and cex
  pl <- par("usr")

  if (is.null(cex)) {
    cex <- par("cex")
  }

  # labels
  if (!is.null(labels)) {
    if (is.numeric(labels)) {
      vals <- round(labels, digits)
    } else {
      vals <- labels
    }
  }

  # bar image
  if (horizontal == TRUE) {
    legend_image <- as.raster(matrix(make_alpha(col, alpha), nrow = 1))
  } else {
    legend_image <- as.raster(matrix(make_alpha(rev(col), alpha), ncol = 1))
  }

  border <- ifelse(is.null(border), NA, border)

  # determining coordinates for elements
  ## bar
  ### bar size
  tw <- diff(pl[1:2])
  th <- diff(pl[3:4])

  w <- tw * width_prop
  h <- th * heigh_prop

  ### bottom left corner
  #### inset adjustment
  insetx <- ifelse(is.null(insetx), inset, insetx)
  insety <- ifelse(is.null(insety), inset, insety)

  #### corner
  if (class(position)[1] == "character") {
    if (position %in% c("topleft", "left", "bottomleft")) {
      xbl <- (pl[1] + (tw * insetx))
    }
    if (position %in% c("bottom", "top", "center")) {
      xbl <- (pl[1] + (tw / 2)) - (w / 2)
    }
    if (position %in% c("bottomright", "topright", "right")) {
      insetx <- ifelse(!horizontal, insetx + 0.08, insetx)
      xbl <- pl[2] - (tw * insetx) - w
    }

    if (position %in% c("left", "right", "center")) {
      ybl <- (pl[3] + (th / 2)) - (h / 2)
    }
    if (position %in% c("bottomright", "bottom", "bottomleft")) {
      insety <- ifelse(!horizontal, insety, insety + 0.03)
      ybl <- pl[3] + (th * insety)
    }
    if (position %in% c("topleft", "top", "topright")) {
      insety <- insety + 0.05
      ybl <- (pl[4] - (th * insety)) - h
    }

    position <- c(xbl, ybl)
  }

  ### bar coordinates
  legend_coord1 <- position[1]
  legend_coord2 <- position[2]
  legend_coord3 <- position[1] + w
  legend_coord4 <- position[2] + h

  xss <- c(position[1], position[1] + w, position[1] + w, position[1])
  yss <- c(position[2], position[2], position[2] + h, position[2] + h)

  ## labels
  if (!is.null(labels)) {
    if (horizontal == FALSE) {
      labels_x <- position[1] + w
      labels_y <- seq(from = position[2], to = (position[2] + h),
                      length.out = length(vals))

      pos <- 4
    } else {
      labels_x <- seq(from = position[1], to = (position[1] + w),
                      length.out = length(vals))
      labels_y <- position[2]

      pos <- 1
    }
  }

  ## title
  if (!is.null(title)) {
    if (horizontal == FALSE) {
      title_coord1 <- position[1]
      title_coord2 <- position[2] + h

      post <- NULL
      adj <- c(0.05, -1)
      toff <- 0.5
    } else {
      title_coord1 <- position[1] + (w / 2)
      title_coord2 <- position[2]

      post <- 3
      adj <- NULL
      toff <- 1
    }
  }


  # plotting bar legend
  ## bar image
  rasterImage(legend_image, legend_coord1, legend_coord2,
              legend_coord3, legend_coord4)
  polygon(x = xss, y = yss, col = NA, border = border)

  ## labels
  if (!is.null(labels)) {
    text(x = labels_x, y = labels_y, labels = vals, pos = pos,
         offset = labels_offset, cex = cex*0.6)
  }

  ## title
  if (!is.null(title)) {
    text(x = title_coord1, y = title_coord2, labels = title,
         cex = cex*0.7, adj = adj, pos = post, offset = toff)
  }
}





# Create a simpler bar legend to be used in plotting functions
bar_legend <- function (value_range, col, alpha = 1, title = NULL, round = 0,
                        label_x = 0.7, labels_y = c(0.2, 0.85),
                        legend_coord = c(0.1, 0.2, 0.3, 0.85),
                        title_coord = c(0.6, 0.525), title_srt = 90,
                        horizontal = FALSE) {
  # Initial tests
  if (missing(value_range)) {
    stop("Argument 'value_range' is required to produce the legend")
  }
  if (missing(col)) {
    stop("Argument 'col' is required to produce the legend")
  }

  # Bar plot
  if (horizontal == TRUE) {
    legend_image <- as.raster(matrix(make_alpha(rev(col), alpha), nrow = 1))
  } else {
    legend_image <- as.raster(matrix(make_alpha(rev(col), alpha), ncol = 1))
  }

  text(x = title_coord[1], y = title_coord[2], labels = title, srt = title_srt)

  if (is.numeric(value_range)) {
    vals <- round(value_range, round)
  } else {
    vals <- value_range
  }

  text(x = label_x, y = labels_y, labels = vals, cex = 0.8)
  rasterImage(legend_image, legend_coord[1], legend_coord[2],
              legend_coord[3], legend_coord[4])
}



# Make colors transparent at distinct levels
make_alpha <- function(col, alpha = 1, names = NULL) {
  rgb_col <- col2rgb(col)
  t_col <- rgb(rgb_col[1, ], rgb_col[2, ], rgb_col[3, ],
               alpha = (alpha * 100) * 255 / 100,
               names = names, maxColorValue = 255)
  return(t_col)
}

