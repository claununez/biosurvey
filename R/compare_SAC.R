
compare_SAC <- function(SAC_selected_sites, element_1, element_2, col_mean1 = "blue",
                        col_CI1 = "lightblue", alpha_mean = 0.7, alpha_CI = 0.2,
                        col_mean2 = "gray15", col_CI2 = "gray65",
                        xlab = "Number of sites", ylab = "Species",
                        line_for_multiple = TRUE, add_legend = TRUE, ...) {
  # Initial tests
  if (missing(SAC_selected_sites)) {
    stop("Argument 'SAC_selected_sites' must be defined.")
  }
  if (missing(element_1)) {
    stop("Argument 'element_1' must be defined.")
  }
  if (missing(element_2)) {
    stop("Argument 'element_2' must be defined.")
  }

  # SACs
  sac1 <- SAC_selected_sites[[element_1]]
  sac2 <- SAC_selected_sites[[element_2]]

  # Names of sites
  sac1nam <- ifelse (is.numeric(element_1), names(SAC_selected_sites)[element_1],
                     element_1)
  sac2nam <- ifelse (is.numeric(element_2), names(SAC_selected_sites)[element_2],
                     element_2)

  sac1nam <- gsub("_", " ", sac1nam)
  sac2nam <- gsub("_", " ", sac2nam)

  # Preparing limits if more than one SAC
  rany1 <- max(unlist(lapply(sac1, function(x) {max(x$richness + (1.96 * x$sd))})))
  rany2 <- max(unlist(lapply(sac2, function(x) {max(x$richness + (1.96 * x$sd))})))

  y_lim <- c(c(0, max(rany1, rany2)))

  # Final colors
  cm1 <- scales::alpha(col_mean1, alpha_mean)
  cm2 <- scales::alpha(col_mean2, alpha_mean)
  ci1 <- scales::alpha(col_CI1, alpha_CI)
  ci2 <- scales::alpha(col_CI2, alpha_CI)

  # Plot
  ## Plot 1
  if (length(sac1) > 1) {
    ### Plot
    if (line_for_multiple == TRUE) {
      pple <- lapply(sac1, function(x) {
        if (x == 1) {
          plot(x, ci.type = "line", ci = 0, col = cm1,
               ylim = y_lim, xlab = xlab, ylab = ylab, ...)
        } else {
          plot(x, ci.type = "line", ci = 0, col = cm1, add = TRUE, ...)
        }
      })
    } else {
      pple <- lapply(sac1, function(x) {
        if (x == 1) {
          plot(x, ci.type = "poly", col =  cm1, ci.lty = 0, ci.col = ci1,
               ylim = y_lim, xlab = xlab, ylab = ylab, ...)
        } else {
          plot(x, ci.type = "poly", col =  cm1, ci.lty = 0, ci.col = ci1,
               add = TRUE, ...)
        }
      })
    }
  } else {
    ## Plot
    plot(sac1[[1]], ci.type = "poly", col =  cm1, ci.lty = 0, ci.col = ci1,
         ylim = y_lim, xlab = xlab, ylab = ylab, ...)
  }

  ## Plot 2
  if (length(sac2) > 1) {
    ### Plot
    if (line_for_multiple == TRUE) {
      pple <- lapply(sac2, function(x) {
        plot(x, ci.type = "line", ci = 0, col = cm2, add = TRUE, ...)
      })
    } else {
      pple <- lapply(sac2, function(x) {
        plot(x, ci.type = "poly", col =  cm2, ci.lty = 0, ci.col = ci2,
             add = TRUE, ...)
      })
    }
  } else {
    ## Plot
    plot(sac2[[1]], ci.type = "poly", col =  cm2, ci.lty = 0, ci.col = ci2,
         add = TRUE, ...)
  }

  if (add_legend == TRUE) {
    legend("bottomright", legend = c(sac1nam, sac2nam), bty = "n", lty = 1,
           col = c(cm1, cm2), cex = 0.8)
  }
}
