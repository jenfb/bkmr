ComputePIPs <- function(fit, rthresh = 10^(-(0:6)), expos.names = NULL, sel = NULL) {
    if(inherits(fit, "stanfit")) {
        rsamps <- rstan::extract(fit)$r
    } else if (inherits(fit, "bkmrfit")) {
        if (is.null(sel)) {
            sel <- with(fit, seq(floor(iter/2) + 1, iter))
        }
        rsamps <- fit$r[sel, ]
    }
    if(is.null(expos.names)) expos.names <- paste0("expos", 1:ncol(rsamps))

    mat <- matrix(NA, length(rthresh), ncol(rsamps), dimnames = list(paste0("rthresh", seq_along(rthresh)), expos.names))
    for(i in seq_along(rthresh)) {
        mat[i, ] <- colMeans(abs(rsamps) >= rthresh[i])
    }
    attr(mat, "rthresh") <- rthresh
    mat
}

#' @import ggplot2
PlotPIPs <- function(pips, rthresh = NULL, plot = TRUE) {
    if (is.null(rthresh)) {
        rthresh <- attr(pips, "rthresh")
        if (is.null(rthresh)) {
            stop("Can't plot unless the thresholds used to compute the PIPs are supplied; use argument rthresh")
        }
    }
    df <- pips %>%
        data.frame %>%
        dplyr::mutate(thresh = rthresh) %>%
        tidyr::gather(exposure, pip, -thresh)

    plt <- ggplot(df, aes(thresh, pip)) +
        geom_path() +
        geom_point() +
        scale_x_log10() +
        coord_cartesian(ylim = c(0, 1)) +
        facet_wrap(~ exposure) +
        labs(x = "Threshold R", y = "Prob(r > R)")

    if(plot) plot(plt)

    ret <- list(df = df, plt = plt)
}
