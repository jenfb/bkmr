#' Compute the posterior mean and variance of \code{h} at a new predictor values
#'
#' Function to approximate the posterior mean and variance as a function of the estimated tau, lambda, beta, and sigsq.eps
#' @inheritParams kmbayes
#' @inheritParams ExtractEsts
#' @param sel A vector selecting which iterations of the BKMR fit should be retained for inference. Currently only implemented for \code{method == "r"}.
#' @export
ComputePostmeanHnew <- function(fit, y = NULL, Z = NULL, X = NULL, Znew, sel = NULL) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
  }
  
    if(is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
    if(class(Znew) == "data.frame") Znew <- data.matrix(Znew)
    if(is.null(dim(X))) X <- matrix(X, ncol=1)

    ests <- ExtractEsts(fit, sel = sel)
    sigsq.eps <- ests$sigsq.eps[, "mean"]
    r <- ests$r[, "mean"]
    beta <- ests$beta[, "mean"]
    lambda <- ests$lambda[, "mean"]

    # if(is.null(data.comps$knots)) {
    n0 <- nrow(Z)
    n1 <- nrow(Znew)
    nall <- n0 + n1
    Kpartall <- makeKpart(r, rbind(Z, Znew))
    Kmat <- exp(-Kpartall)
    Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
    Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
    Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]

    Kpart <- makeKpart(r, Z)
    V <- diag(1, nrow(Z), nrow(Z)) + lambda[1]*exp(-Kpart)
    cholV <- chol(V)
    Vinv <- chol2inv(cholV)

    lamK10Vinv <- lambda[1]*Kmat10 %*% Vinv
    Sigma.hnew <- lambda[1]*sigsq.eps*(Kmat1 - lamK10Vinv %*% t(Kmat10))
    mu.hnew <- lamK10Vinv %*% (y - X%*%beta)
    # } else {
    # stop("GPP not yet implemented")
    # }

    ret <- list(postmean = drop(mu.hnew), postvar = Sigma.hnew)
    ret
}
