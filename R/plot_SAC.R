


plot_SAC <- function(SAC_selected_sites, col_mean = "blue", col_CI = "lightblue",
                      alpha_mean = 0.7, alpha_CI = 0.2, xlab = "Number of sites",
                      ylab = "Species", ...) {

  # Inital tests
  if (missing(PAM_selection)) {
    stop("Argument 'PAM_selection' must be defined.")
  }
  if (missing(SAC1)) {
    stop("Argument 'SAC1' must be defined.")
  }
  if (missing(SAC2)) {
    stop("Argument 'SAC2' must be defined.")
  }

  # Needed library
  suppressPackageStartupMessages(library(vegan))

  # Plotting
  if (length(SAC_selected_sites) == 1) {
    if (length(SAC_selected_sites[[1]]) > 1) {
      ## Preparing limits if more than one SAC
      maxy <- max(sapply(SAC_selected_sites[[1]], function(x) {max(x$richness)}))
      y_lim <- c(0, maxy)

      ## Plot
      pple <- lapply(SAC_selected_sites[[1]], function(x) {
        if (x == 1) {
          plot(x, ci.type = "line", ci = 0,
               col = scales::alpha(col_mean, alpha_mean),
               ylim = y_lim, xlab = xlab, ylab = ylab, ...)
        } else {
          plot(x, ci.type = "line", ci = 0,
               col = scales::alpha(col_mean, alpha_mean), add = TRUE, ...)
        }
      })
    } else {
      ## Limits for one SAC
      y_lim <- c(0, max(SAC_selected_sites[[1]][[1]]$richness))

      ## Plot
      plot(SAC_selected_sites[[1]][[1]], ci.type = "poly",
           col =  scales::alpha(col_mean, alpha_mean),
           ci.lty = 0, ci.col = scales::alpha(col_CI, alpha_CI),
           ylim = y_lim, xlab = xlab, ylab = ylab, ...)
    }

  } else {
    ## Defining par settings
    nl <- length(y)
    if (missing(nc)) {
      nc <- ceiling(sqrt(nl))
    } else {
      nc <- max(1, min(nl, round(nc)))
    }
    if (missing(nr)) {
      nr <- ceiling(nl / nc)
    } else {
      nr <- max(1, min(nl, round(nr)))
      nc <- ceiling(nl / nr)
    }

    par(mfrow = c(nr, nc), cex = par_cex)

    ## Plots in loop
    for (i in 1:length(SAC_selected_sites)) {
      if (length(SAC_selected_sites[[i]]) > 1) {
        ## Preparing limits if more than one SAC
        maxy <- max(sapply(SAC_selected_sites[[i]], function(x) {max(x$richness)}))
        y_lim <- c(0, maxy)

        ## Plot
        pple <- lapply(SAC_selected_sites[[i]], function(x) {
          if (x == 1) {
            plot(x, ci.type = "line", ci = 0,
                 col = scales::alpha(col_mean, alpha_mean),
                 ylim = y_lim, xlab = xlab, ylab = ylab, ...)
          } else {
            plot(x, ci.type = "line", ci = 0,
                 col = scales::alpha(col_mean, alpha_mean), add = TRUE, ...)
          }
        })
      } else {
        ## Limits for one SAC
        y_lim <- c(0, max(SAC_selected_sites[[i]][[1]]$richness))

        ## Plot
        plot(SAC_selected_sites[[i]][[1]], ci.type = "poly",
             col =  scales::alpha(col_mean, alpha_mean),
             ci.lty = 0, ci.col = scales::alpha(col_CI, alpha_CI),
             ylim = y_lim, xlab = xlab, ylab = ylab, ...)
      }
    }
  }
}
