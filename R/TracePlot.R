#' Trace plot
#'
#' @param fit 
#' @param par 
#' @param comp 
#' @param sel 
#' @param main 
#' @param xlab 
#' @param ylab 
#' @param ... 
#'
#' @export
TracePlot <- function(fit, par, comp = 1, sel = NULL, main = "", xlab = "iteration", ylab = "parameter value", ...) {
    samps <- ExtractSamps(fit, sel = sel)[[par]]
    if (!is.null(ncol(samps))) {
        nm <- colnames(samps)[comp]
        samps <- samps[, comp]
    } else {
        nm <- par
    }
    main <- paste0(main, "\n(", nm, " = ", format(mean(samps), digits = 2), ")")
    plot(samps, type = "l", main = main, xlab = xlab,  ylab = ylab, ...)
    abline(h = mean(samps), col = "blue", lwd = 2)
}
