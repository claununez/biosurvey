


plot_SACs <- function(PAM_selection, SAC1, SAC2, main = NULL, col_mean1 = "gray45", 
                      col_mean2 = "blue", col_CI1 = "gray65", col_CI2 = "lightblue", 
                      alpha_mean = 0.6, alpha_CI = 0.3, xlab = "Number of sites", 
                      ylab = "Species", ...) {
  
  if (missing(PAM_selection)) {
    stop("Argument 'PAM_selection' must be defined.")
  }
  if (missing(SAC1)) {
    stop("Argument 'SAC1' must be defined.")
  }
  if (missing(SAC2)) {
    stop("Argument 'SAC2' must be defined.")
  }
  if (is.null(main)) {
    main <- ""
  }
  SACl <- list(SAC1, SAC2)
  cols <- c(col1, col2)
  colsl <- c(col_lim1, col_lim2)
  
  if (!is.null(is_list)) {
    if (is_list == 1) {
      or <- 1:2
    } else {
      or <- 2:1
    }
    
    SACl <- SACl[or]
    cols <- cols[or]
    colsl <- colsl[or]
    
    if (is.null(y_lim)) {
      all_max <- c(max(SACl[[2]]$richness), 
                   max(sapply(SACl[[1]], function(x) {max(x$richness)})))
      y_lim <- c(0, max(all_max))
    }
  } else {
    if (is.null(y_lim)) {y_lim <- c(0, max(sapply(SACl, function(x) {max(x$richness)})))}
  }
  
  par(mar = par_mar, cex = par_cex)
  if (!is.null(is_list)) {
    pple <- lapply(1:length(SACl[[1]]), function(x) {
      if (x == 1) {
        plot(SACl[[1]][[x]], ci.type = "line", ci = 0, 
             col = scales::alpha(cols[1], alpha_mean), 
             ylim = y_lim, xlab = xlab, ylab = ylab, main = main, ...)
      } else {
        plot(SACl[[1]][[x]], ci.type = "line", ci = 0, 
             col = scales::alpha(cols[1], alpha_mean), add = TRUE, ...)
      }
    })
  } else {
    plot(SACl[[1]], ci.type = "poly", col =  scales::alpha(cols[1], alpha_mean), 
         ci.lty = 0, ci.col = scales::alpha(colsl[1], alpha_CI), 
         ylim = y_lim, xlab = xlab, ylab = ylab, main = main, ...)  
  }
  plot(SACl[[2]], ci.type = "poly", col =  scales::alpha(cols[2], alpha_mean), 
       ci.lty = 0, ci.col = scales::alpha(colsl[2], alpha_CI), add = TRUE, ...)
}