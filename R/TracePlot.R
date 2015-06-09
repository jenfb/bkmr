#' Trace plot
#'
#' @export
TracePlot <- function(fit, par, comp = 1, main = "", ...) {
    if(inherits(fit, "stanfit")) {
        warning("you can use rstan::traceplot() instead")
    }
    samps <- ExtractSamps(fit)[[par]]
    if (!is.null(ncol(samps))) {
        nm <- colnames(samps)[comp]
        samps <- samps[, comp]
    } else {
        nm <- par
    }
    main <- paste0(main, "\n(", nm, " = ", format(mean(samps), digits = 2), ")")
    plot(samps, type = "l", ylab = nm, main = main, ...)
    abline(h = mean(samps), col = "blue", lwd = 2)
}
