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
#' library(rstan)
#' set.seed(1)
#' dat <- SimData(n = 50, M = 4, sigsq.true = 0.5, beta.true = 2)
#' y <- dat$y
#' Z <- dat$Z
#' X <- dat$X
#' set.seed(111)
#' runtime <- system.time({
#'     fitstanex1 <- StanBKMR(y = y, Z = Z, X = X, iter = 400, chains = 3)
#' })
#' print(fitstanex1, par = c("beta", "sigma_sq", "h[1]", "h[100]", "lp__"))
#' }
StanBKMR <- function(y, Z, X, file = NULL, ...) {
    require(rstan)

    if (is.null(file)) {
        file <- system.file("stan", "gp-multi-fit.stan", package = "bkmr")
    }

    data <- list(D = ncol(Z), K = ncol(X), N = length(y), Z = Z, X = X, y = y)

    fitstan <- rstan::stan(file = file, data = data, ...)
}




