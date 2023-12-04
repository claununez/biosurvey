#' Plotting lists of species accumulation curves
#'
#' @description Creates species accumulation curve plots (one or multiple
#' panels) from information contained in lists obtained with the function
#' \code{\link{selected_sites_SAC}}.
#'
#' @param SAC_selected_sites nested list of "\code{specaccum}" objects obtained
#' with function \code{\link{selected_sites_SAC}}.
#' @param col_mean (character) color for mean value of curve. Default = "blue".
#' @param col_CI (character) color for confidence interval region for the curve.
#' Default = "lightblue".
#' @param alpha_mean (numeric) alpha level for line representing the mean,
#' values from 0 to 1. Default = 0.7. Values close to 0 increase transparency.
#' @param alpha_CI (numeric) alpha level for the region representing the
#' confidence interval. Default = 0.2.
#' @param xlab (character) label for x-axis of plot. Default = "Number of
#' sites".
#' @param ylab (character) label for y-axis of plot. Default = "Species".
#' @param line_for_multiple (logical) whether to plot SACs only as lines when
#' multiple objects are in one or more of the internal lists in
#' \code{SAC_selected_sites}. Default = TRUE.
#' @param main (character) title or titles for plots. The default, NULL, adds
#' titles according to names of elements in \code{SAC_selected_sites}.
#' @param ... other arguments to be passed to plot method for objects of class
#' "\code{specaccum}".
#'
#' @return
#' A plot of "\code{specaccum}" objects. Multiple panels will be plotted
#' if \code{SAC_selected_sites} list contains more than one element.
#'
#' @usage
#' plot_SAC(SAC_selected_sites, col_mean = "blue", col_CI = "lightblue",
#'          alpha_mean = 0.7, alpha_CI = 0.2, xlab = "Number of sites",
#'          ylab = "Species", line_for_multiple = TRUE, main = NULL, ...)
#'
#' @export
#' @import vegan
#'
#' @examples
#' # Data
#' b_pam <- read_PAM(system.file("extdata/b_pam.rds",
#'                               package = "biosurvey"))
#' m_selection <- read_master(system.file("extdata/m_selection.rds",
#'                                        package = "biosurvey"))
#'
#' # Subsetting base PAM according to selections
#' sub_pam_all <- subset_PAM(b_pam, m_selection, selection_type = "all")
#'
#' SACs <- selected_sites_SAC(PAM_subset = sub_pam_all, selection_type = "all")
#'
#' # Plotting
#' plot_SAC(SACs)

plot_SAC <- function(SAC_selected_sites, col_mean = "blue",
                     col_CI = "lightblue", alpha_mean = 0.7, alpha_CI = 0.2,
                     xlab = "Number of sites", ylab = "Species",
                     line_for_multiple = TRUE, main = NULL, ...) {

  # Initial tests
  if (missing(SAC_selected_sites)) {
    stop("Argument 'SAC_selected_sites' must be defined.")
  }

  lsac <- length(SAC_selected_sites)

  # Par settings
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  # Final colors
  cm <- make_alpha(col_mean, alpha_mean)
  cci <- make_alpha(col_CI, alpha_CI)

  # Y limits
  maxy <- max(unlist(lapply(SAC_selected_sites, function(w) {
    sapply(w, function(x) {max(x$richness + (1.96 * x$sd))})
  })))

  y_lim <- c(0, maxy)


  # Mains
  if (is.null(main)) {
    mains <- gsub("_", " ", names(SAC_selected_sites))
  } else {
    ## Options if not null
    if (class(main)[1] != "character") {
      main <- as.character(main)
    }

    if (length(main) == 1) {
      mains <- rep(main, lsac)
    } else {
      if (length(main) != lsac) {
        message("length of 'main' does not coincide with length of 'SAC_selected_sites'",
                "\n'main' elements in 'main' will be recycled.")
        mains <- rep(main, lsac)[1:lsac]
      } else {
        mains <- main
      }
    }
  }

  # Plotting
  if (lsac == 1) {
    if (length(SAC_selected_sites[[1]]) > 1) {
      ## Plot
      if (line_for_multiple == TRUE) {
        pple <- lapply(1:length(SAC_selected_sites[[1]]), function(x) {
          if (x == 1) {
            plot(SAC_selected_sites[[1]][[x]], ci.type = "line", ci = 0,
                 col = cm, ylim = y_lim, xlab = xlab, ylab = ylab,
                 main = mains[1], ...)
          } else {
            plot(SAC_selected_sites[[1]][[x]], ci.type = "line", ci = 0,
                 col = cm, add = TRUE, ...)
          }
        })
      } else {
        pple <- lapply(1:length(SAC_selected_sites[[1]]), function(x) {
          if (x == 1) {
            plot(SAC_selected_sites[[1]][[x]], ci.type = "poly", col =  cm,
                 ci.lty = 0, ci.col = cci, ylim = y_lim, xlab = xlab,
                 ylab = ylab, main = mains[1], ...)
          } else {
            plot(SAC_selected_sites[[1]][[x]], ci.type = "poly", col =  cm,
                 ci.lty = 0, ci.col = cci,
                 add = TRUE, ...)
          }
        })
      }
    } else {
      ## Plot
      plot(SAC_selected_sites[[1]][[1]], ci.type = "poly", col =  cm,
           ci.lty = 0, ci.col = cci, ylim = y_lim, xlab = xlab, ylab = ylab,
           main = mains[1], ...)
    }

  } else {
    ## Defining new par settings
    nl <- length(SAC_selected_sites)
    nc <- ceiling(sqrt(nl))
    nr <- ceiling(nl / nc)

    par(mfrow = c(nr, nc))

    ## Plots in loop
    for (i in 1:nl) {
      if (length(SAC_selected_sites[[i]]) > 1) {
        ## Plot
        if (line_for_multiple == TRUE) {
          pple <- lapply(1:length(SAC_selected_sites[[i]]), function(x) {
            if (x == 1) {
              plot(SAC_selected_sites[[i]][[x]], ci.type = "line", ci = 0,
                   col = cm, ylim = y_lim, xlab = xlab, ylab = ylab,
                   main = mains[i], ...)
            } else {
              plot(SAC_selected_sites[[i]][[x]], ci.type = "line", ci = 0,
                   col = cm, add = TRUE, ...)
            }
          })
        } else {
          pple <- lapply(1:length(SAC_selected_sites[[i]]), function(x) {
            if (x == 1) {
              plot(SAC_selected_sites[[i]][[x]], ci.type = "poly", col =  cm,
                   ci.lty = 0, ci.col = cci, ylim = y_lim, xlab = xlab,
                   ylab = ylab, main = mains[i], ...)
            } else {
              plot(SAC_selected_sites[[i]][[x]], ci.type = "poly", col =  cm,
                   ci.lty = 0, ci.col = cci, add = TRUE, ...)
            }
          })
        }
      } else {
        ## Plot
        plot(SAC_selected_sites[[i]][[1]], ci.type = "poly", col =  cm,
             ci.lty = 0, ci.col = cci, ylim = y_lim, xlab = xlab, ylab = ylab,
             main = mains[i], ...)
      }
    }
  }
}
