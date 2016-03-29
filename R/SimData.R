#' Simulate dataset
#'
#' Simulate predictor, covariate, and outcome data
#'
#' @export
#'
#' @param n Number of observations
#' @param M Number of predictor variables to generate
#' @param sigsq.true Variance of normally distributed residual error
#' @param beta.true Coefficient on the covariate
#' @param hfun An integer from 1 to 3 identifying which predictor-response function to generate
#' @examples
#' set.seed(5)
#' dat <- SimData()
#' @details 
#' \itemize{
#'  \item{"hfun = 1"}{A nonlinear function of the first predictor}
#'  \item{"hfun = 2"}{A linear function of the first two predictors and their product term}
#'  \item{"hfun = 3"}{A nonlinear and nonadditive function of the first two predictor variables}
#' }
SimData <- function(n = 100, M = 5, sigsq.true = 0.5,
                     beta.true = 2, hfun = 3) {
  
  stopifnot(n > 0, M > 0, sigsq.true >= 0)
  
  if (hfun == 1) {
    HFun <- function(z, ind = 1) 4*plogis(z[ind], 0, 0.3)
  } else if (hfun == 2) {
    HFun <- function(z, ind1 = 1, ind2 = 2) 1/4*(z[ind1] + z[ind2] + 1/2*z[ind1]*z[ind2])
  } else if (hfun == 3) {
    HFun <- function(z, ind1 = 1, ind2 = 2) 4*plogis(1/4*(z[ind1] + z[ind2] + 1/2*z[ind1]*z[ind2]), 0, 0.3)
  } else {
    stop("hfun must be an integer from 1 to 3")
  }
  
  Z <- matrix(rnorm(n * M), n, M,
                  dimnames = list(NULL, paste0("z", 1:M)))
  X <- cbind(3*cos(Z[, 1]) + 2*rnorm(n))
  eps <- rnorm(n, sd = sqrt(sigsq.true))
  h <- apply(Z, 1, HFun)
  y <- drop(X * beta.true + h + eps)
  
  dat <- list(n = n, M = M, sigsq.true = sigsq.true, beta.true = beta.true, Z = Z, h = h, X = X, y = y, hfun = hfun, HFun = HFun)
}
