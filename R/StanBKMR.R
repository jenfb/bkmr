#' Wrapper function to fit BKMR using Stan
#'
#' @export
StanBKMR <- function(y, expos, covar, file, ...) {
    if (missing(file)) {
        file <- system.file("stan", "gp-multi-fit.stan", package = "bkmr")
    }

    data <- list(D = ncol(expos), K = ncol(covar), N = length(y), expos = expos, covar = covar, y = y)

    fitstan <- stan(file = system.file("stan", "gp-multi-fit.stan", package = "bkmr"), data = data, ...)
}




