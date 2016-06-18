HFun1 <- function(z, ind = 1) 4*plogis(z[ind[1]], 0, 0.3)
HFun2 <- function(z, ind = 1:2) 1/4*(z[ind[1]] + z[ind[2]] + 1/2*z[ind[1]]*z[ind[2]])
HFun3 <- function(z, ind = 1:2) 4*plogis(1/4*(z[ind[1]] + z[ind[2]] + 1/2*z[ind[1]]*z[ind[2]]), 0, 0.3)

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
#' @param Zgen Method for generating the matrix Z of exposure variables, taking one of the values c("unif", "corr", "realistic")
#' @param ind select which predictor(s) will be included in the \code{h} function; how many predictors that can be included will depend on which \code{h} function is being used.
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
                    beta.true = 2, hfun = 3, Zgen = "unif", ind = 1:2) {
  
  stopifnot(n > 0, M > 0, sigsq.true >= 0)
  
  if (hfun == 1) {
    HFun <- HFun1
  } else if (hfun == 2) {
    HFun <- HFun2
  } else if (hfun == 3) {
    HFun <- HFun3
  } else {
    stop("hfun must be an integer from 1 to 3")
  }
  
  if (Zgen == "unif") {
    Z <- matrix(rnorm(n * M), n, M)
  } else if (Zgen == "corr") {
    if (M < 3) {
      stop("M must be an integer > 2 for Zgen = 'corr'")
    }
    Sigma <- diag(1, M, M)
    Sigma[1,3] <- Sigma[3,1] <- 0.95
    Sigma[2,3] <- Sigma[3,2] <- 0.3
    Sigma[1,2] <- Sigma[2,1] <- 0.1
    Z <- MASS::mvrnorm(n = n, mu = rep(0, M), Sigma = Sigma)
  } else if (Zgen == "realistic") {
    VarRealistic <- structure(c(0.72, 0.65, 0.45, 0.48, 0.08, 0.14, 0.16, 0.42, 0.2, 
                                0.11, 0.35, 0.1, 0.11, 0.65, 0.78, 0.48, 0.55, 0.06, 0.09, 0.17, 
                                0.2, 0.16, 0.11, 0.32, 0.12, 0.12, 0.45, 0.48, 0.56, 0.43, 0.11, 
                                0.15, 0.23, 0.25, 0.28, 0.16, 0.31, 0.15, 0.14, 0.48, 0.55, 0.43, 
                                0.71, 0.2, 0.23, 0.32, 0.22, 0.29, 0.14, 0.3, 0.22, 0.18, 0.08, 
                                0.06, 0.11, 0.2, 0.95, 0.7, 0.45, 0.22, 0.29, 0.16, 0.24, 0.2, 
                                0.13, 0.14, 0.09, 0.15, 0.23, 0.7, 0.8, 0.36, 0.3, 0.35, 0.13, 
                                0.23, 0.17, 0.1, 0.16, 0.17, 0.23, 0.32, 0.45, 0.36, 0.83, 0.24, 
                                0.37, 0.2, 0.36, 0.34, 0.25, 0.42, 0.2, 0.25, 0.22, 0.22, 0.3, 
                                0.24, 1.03, 0.41, 0.13, 0.39, 0.1, 0.1, 0.2, 0.16, 0.28, 0.29, 
                                0.29, 0.35, 0.37, 0.41, 0.65, 0.18, 0.3, 0.18, 0.16, 0.11, 0.11, 
                                0.16, 0.14, 0.16, 0.13, 0.2, 0.13, 0.18, 0.6, 0.18, 0.13, 0.08, 
                                0.35, 0.32, 0.31, 0.3, 0.24, 0.23, 0.36, 0.39, 0.3, 0.18, 0.79, 
                                0.42, 0.12, 0.1, 0.12, 0.15, 0.22, 0.2, 0.17, 0.34, 0.1, 0.18, 
                                0.13, 0.42, 1.27, 0.1, 0.11, 0.12, 0.14, 0.18, 0.13, 0.1, 0.25, 
                                0.1, 0.16, 0.08, 0.12, 0.1, 0.67), .Dim = c(13L, 13L))
    if (M > ncol(VarRealistic)) {
      stop("Currently can only generate exposure data based on a realistic correlation structure with M = 13 or fewer. Please set M = 13 or use Zgen = 'unif'")
    } else if (M <= 13) {
      Sigma <- VarRealistic[1:M, 1:M]
    }
    Z <- MASS::mvrnorm(n = n, mu = rep(0, M), Sigma = Sigma)
  }
  colnames(Z) <- paste0("z", 1:M)
  
  X <- cbind(3*cos(Z[, 1]) + 2*rnorm(n))
  eps <- rnorm(n, sd = sqrt(sigsq.true))
  h <- apply(Z, 1, HFun)
  y <- drop(X * beta.true + h + eps)
  
  dat <- list(n = n, M = M, sigsq.true = sigsq.true, beta.true = beta.true, Z = Z, h = h, X = X, y = y, hfun = hfun, HFun = HFun)
}
