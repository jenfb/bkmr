GenCovarDefault <- function(z) {
    3*cos(z[1]) + 2*rnorm(1)
}

GenExposDefault <- function(n, M) {
    matrix(rnorm(n * M), n, M,
           dimnames = list(NULL, paste0("expos", 1:M)))
}

HFun1 <- function(z, ind = 1) 4*plogis(z[ind], 0, 0.3)

#' Simulate dataset
#'
#' Simulate dataset
#'
#' @export
#'
#' @param GenExpos A function that takes as input a vector of length \code{M}
#' @param GenCovar A function that takes as input a vector of length \code{M}
#' @param HFun A function that takes as input a vector of length \code{M} and outputs a scalar
SimData1 <- function(n = 100, M = 5, sigsq.true = 0.5,
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
