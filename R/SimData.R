GenCovarDefault <- function(z) {
    3*cos(z[1]) + 2*rnorm(1)
}

GenExposDefault <- function(n, M) {
    matrix(rnorm(n * M), n, M,
           dimnames = list(NULL, paste0("expos", 1:M)))
}

HFun1 <- function(z, ind = 1) 4*plogis(z[ind], 0, 0.3)

HFun2 <- function(z, ind1 = 1, ind2 = 2) 1/4*(z[ind1] + z[ind2] + 1/2*z[ind1]*z[ind2])

HFun3 <- function(z, ind1 = 1, ind2 = 2) 4*plogis(1/4*(z[ind1] + z[ind2] + 1/2*z[ind1]*z[ind2]), 0, 0.3)

#' Simulate dataset
#'
#' Simulate dataset
#'
#' @export
#'
#' @param GenExpos A function that takes as input a vector of length \code{M} Defaults to \code{\link{GenExposDefault}}.
#' @param GenCovar A function that takes as input a vector of length \code{M}. Defaults to \code{\link{GenCovarDefault}}.
#' @param n Number of observations
#' @param M Number of exposure variables to generate
#' @param sigsq.true Variance of normally distributed residual error
#' @param beta.true Coefficient on the covariate
#' @param HFun A function that takes as input a vector of length \code{M} and outputs a scalar
#' @examples
#' set.seed(5)
#' dat <- SimData()
SimData <- function(n = 100, M = 5, sigsq.true = 0.5,
                     beta.true = 2, GenCovar = NULL,
                     GenExpos = NULL, HFun = NULL) {
    if (is.null(GenCovar)) {
        GenCovar <- GenCovarDefault
    }
    if (is.null(GenExpos)) {
        GenExpos <- GenExposDefault
    }
    if (is.null(HFun)) {
        HFun <- HFun1
    }
    expos <- GenExpos(n = n, M = M)
    eps <- rnorm(n, sd = sqrt(sigsq.true))
    h <- apply(expos, 1, HFun)
    covar <- cbind(apply(expos, 1, GenCovar))
    if (length(beta.true) == 1) {
        beta.true <- rep(beta.true, ncol(covar))
    }
    y <- drop(covar %*% beta.true + h + eps)

    dat <- list(n = n, M = M, sigsq.true = sigsq.true, beta.true = beta.true, expos = expos, h = h, covar = covar, y = y)
}
