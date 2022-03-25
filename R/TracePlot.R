#' Trace plot
#'
#' @inheritParams ExtractEsts
#' @param par which parameter to plot
#' @param comp which component of the parameter vector to plot
#' @param main title
#' @param xlab x axis label
#' @param ylab y axis label
#' @param ... other arguments to pass onto the plotting function
#' @export
#' @import graphics
#' @details For guided examples, go to \url{https://jenfb.github.io/bkmr/overview.html}
#' 
#' @return No return value, generates plot
#' 
#' @examples
#' ## First generate dataset
#' set.seed(111)
#' dat <- SimData(n = 50, M = 4)
#' y <- dat$y
#' Z <- dat$Z
#' X <- dat$X
#' 
#' ## Fit model with component-wise variable selection
#' ## Using only 100 iterations to make example run quickly
#' ## Typically should use a large number of iterations for inference
#' set.seed(111)
#' fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE)
#' 
#' TracePlot(fit = fitkm, par = "beta")
#' TracePlot(fit = fitkm, par = "sigsq.eps")
#' TracePlot(fit = fitkm, par = "r", comp = 1)
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
