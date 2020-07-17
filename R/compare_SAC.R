#' Comparative plots of species accumulation curves
#'
#' @description Creates comparative plots of two species accumulation curves
#' from information contained in lists obtained with the function
#' \code{\link{selected_sites_SAC}}.
#'
#' @param SAC_selected_sites nested list of "\code{specaccum}" objects obtained
#' with function \code{\link{selected_sites_SAC}}.
#' @param element_1 (numeric or character) index of position or name of the first
#' element (type of selection) present in \code{SAC_selected_sites} to be plotted.
#' @param element_2 (numeric or character) index of position or name of the second
#' element (type of selection) present in \code{SAC_selected_sites} to be plotted.
#' @param col_mean1 (character) color for mean value of curve in \code{element_1};
#' default = "blue".
#' @param col_mean2 (character) color for mean value of curve in \code{element_2};
#' default = "gray15".
#' @param col_CI1 (character) color for confidence interval region for the curve
#' in \code{element_1}; default = "lightblue".
#' @param col_CI2 (character) color for confidence interval region for the curve
#' in \code{element_2}; default = "gray65".
#' @param lty1 type of line for \code{element_1}. See lty in \code{\link{par}}.
#' @param lty2 type of line for \code{element_2}.
#' @param alpha_mean (numeric) alpha level for line representing the mean, values
#' from 0 to 1; default = 0.9. Values close to 0 increase transparency.
#' @param alpha_CI (numeric) alpha level for the region representing the confidence
#' interval; default = 0.3.
#' @param xlab (character) label for x axis of plot; default = "Number of sites".
#' @param ylab (character) label for y axis of plot; default = "Species".
#' @param line_for_multiple (logical) whether to plot SACs only as lines when
#' multiple objects are in one or more of the internal lists in
#' \code{SAC_selected_sites}. Default = TRUE.
#' @param add_legend (logical) whether to add default legend to plot; default =
#' TRUE.
#' @param ... other arguments to be passed to plot method for objects of class
#' "\code{specaccum}".
#'
#' @return
#' A comparative plot of two species "\code{specaccum}" objects done based on what
#' is defined in \code{element_1} and  \code{element_2}.
#'
#' @usage
#' compare_SAC(SAC_selected_sites, element_1, element_2, col_mean1 = "blue",
#'             col_CI1 = "lightblue", col_mean2 = "gray15", col_CI2 = "gray65",
#'             alpha_mean = 0.7, alpha_CI = 0.2, xlab = "Number of sites",
#'             ylab = "Species", line_for_multiple = TRUE, add_legend = TRUE, ...)
#'
#' @export
#' @importFrom scales alpha
#' @importFrom graphics legend
#' @import vegan
#'
#' @examples
#' # Data
#' data("b_pam", package = "biosurvey")
#' data("m_selection", package = "biosurvey")
#'
#' # Subsetting base PAM according to selections
#' sub_pam_all <- subset_PAM(b_pam, m_selection, selection_type = "all")
#'
#' SACs <- selected_sites_SAC(PAM_subset = sub_pam_all, selection_type = "all")
#'
#' compare_SAC(SAC_selected_sites = SACs, element_1 = 1, element_2 = 2)

compare_SAC <- function(SAC_selected_sites, element_1, element_2, col_mean1 = "blue",
                        col_CI1 = "lightblue", col_mean2 = "gray15",
                        col_CI2 = "gray65", lty1 = 1, alpha_mean = 0.9, alpha_CI = 0.3,
                        lty2 = 2, xlab = "Number of sites", ylab = "Species",
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
          plot(x, ci.type = "line", ci = 0, col = cm1, lty = lty1,
               ylim = y_lim, xlab = xlab, ylab = ylab, ...)
        } else {
          plot(x, ci.type = "line", ci = 0, col = cm1, lty = lty1, add = TRUE, ...)
        }
      })
    } else {
      pple <- lapply(sac1, function(x) {
        if (x == 1) {
          plot(x, ci.type = "poly", col =  cm1, ci.lty = 0, ci.col = ci1,
               lty = lty1, ylim = y_lim, xlab = xlab, ylab = ylab, ...)
        } else {
          plot(x, ci.type = "poly", col =  cm1, ci.lty = 0, ci.col = ci1,
               lty = lty1, add = TRUE, ...)
        }
      })
    }
  } else {
    ## Plot
    plot(sac1[[1]], ci.type = "poly", col =  cm1, ci.lty = 0, ci.col = ci1,
         lty = lty1, ylim = y_lim, xlab = xlab, ylab = ylab, ...)
  }

  ## Plot 2
  if (length(sac2) > 1) {
    ### Plot
    if (line_for_multiple == TRUE) {
      pple <- lapply(sac2, function(x) {
        plot(x, ci.type = "line", ci = 0, col = cm2, lty = lty2, add = TRUE, ...)
      })
    } else {
      pple <- lapply(sac2, function(x) {
        plot(x, ci.type = "poly", col =  cm2, ci.lty = 0, ci.col = ci2,
             lty = lty2, add = TRUE, ...)
      })
    }
  } else {
    ## Plot
    plot(sac2[[1]], ci.type = "poly", col =  cm2, ci.lty = 0, ci.col = ci2,
         lty = lty2, add = TRUE, ...)
  }

  if (add_legend == TRUE) {
    legend("bottomright", legend = c(sac1nam, sac2nam), bty = "n",
           lty = c(lty1, lty2), col = c(cm1, cm2), cex = 0.8)
  }
}
