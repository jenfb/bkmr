#' Wrapper function to fit BKMR using Stan
#'
#' Fits Bayesian kernel machine regression (BKMR) for cross-sectional data using \code{\link[rstan]{stan}}.
#'
#' @inheritParams kmbayes
#' @param Optional path to the file containing Stan model code.
#' @param ...
#'
#' @return Object of class \code{"\link[=stanfit-class]{stanfit}"}
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' dat <- SimData1(n = 50, M = 4, sigsq.true = 0.5, beta.true = 2)
#' y <- dat$y
#' expos <- dat$expos
#' covar <- dat$covar
#' set.seed(111)
#' runtime <- system.time({
#'     fitstanex1 <- StanBKMR(y = y, expos = expos, covar = covar, iter = 400, chains = 3)
#' })
#' print(fitstanex1, par = c("beta", "sigma_sq", "h[1]", "h[100]", "lp__"))
#' }
StanBKMR <- function(y, expos, covar, file = NULL, ...) {
    require(rstan)

    if (is.null(file)) {
        file <- system.file("stan", "gp-multi-fit.stan", package = "bkmr")
    }

    data <- list(D = ncol(expos), K = ncol(covar), N = length(y), expos = expos, covar = covar, y = y)

    fitstan <- rstan::stan(file = file, data = data, ...)
}




