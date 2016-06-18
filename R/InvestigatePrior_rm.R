#' Investigate prior
#'
#' Investigate the impact of the \code{r[m]} parameters on the smoothness of the exposure-response function \code{h(z[m])}.
#'
#' @inheritParams kmbayes
#' @param ngrid Number of grid points over which to plot the exposure-response function
#' @param q.seq Sequence of values corresponding to different degrees of smoothness in the estimated exposure-response function. A value of q corresponds to fractions of the range of the data over which there is a decay in the correlation \code{cor(h[i],h[j])} between two subjects by 50\code{\%}.
#' @param r.seq sequence of values at which to fix \code{r} for estimating the exposure-response function
#' @param verbose TRUE or FALSE: flag indicating whether to print to the screen which exposure variable and q value has been completed
#' @param Drange the range of the \code{z_m} data over which to apply the values of \code{q.seq}. If not specified, will be calculated as the maximum of the ranges of \code{z_1} through \code{z_M}.
#' @export
#' @import nlme
InvestigatePrior <- function(y, Z, X, ngrid = 50, q.seq = c(2, 1, 1/2, 1/4, 1/8, 1/16), r.seq = NULL, Drange = NULL, verbose = FALSE) {
  
  if (is.null(r.seq)) {
    if (is.null(Drange)) {
      zranges <- diff(apply(Z, 2, range))
      Drange <- max(zranges)
    }
    r.seq <- -log(1 - 0.50)/(q.seq * Drange)^2
  } else {
    if (!is.null(q.seq)) {
      warning("Both 'q.seq' and 'r.seq' are specified; values of 'q.seq' will be ignored.")
    }
  }
  
  Znew.mat <- matrix(NA, ngrid, ncol(Z))
  preds <- vector("list", ncol(Z))
  resids <- vector("list", ncol(Z))
  h.hat.ests <- vector("list", ncol(Z))
  for(i in 1:ncol(Z)) {
    preds[[i]] <- matrix(NA, ngrid, length(q.seq))
    resids[[i]] <- matrix(NA, nrow(Z), length(q.seq))
    h.hat.ests[[i]] <- matrix(NA, nrow(Z), length(q.seq))
  }
  for(i in 1:ncol(Z)) {
    Zi <- cbind(Z[, i])
    Znew <- cbind(seq(min(Zi), max(Zi), length.out = ngrid))
    Znew.mat[, i] <- Znew
    n0 <- nrow(Zi)
    In <- diag(1,n0,n0)
    n1 <- nrow(Znew)
    Zall <- rbind(Zi, Znew)
    nall <- n0+n1
    for(j in seq_along(r.seq)) {
      r <- r.seq[j]
      Kpart <- as.matrix(stats::dist(Zall))^2
      Kmat <- exp(-r*Kpart)
      K <- Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
      Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
      Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]
      U <- try(t(chol(K)), silent=TRUE)
      # all.equal(K, U %*% t(U))
      if(class(U) == "try-error") {
        sigsvd <- svd(K)
        U <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
        # all.equal(K, U %*% t(U), check.attributes=FALSE)
      }
      
      group <- rep(1, n0)
      fit <- lme(y ~ -1+X, random = list(group = pdIdent(~ -1+U)))
      #data.frame(sig=(fit$sigma)^2, tau=as.numeric(VarCorr(fit)[1,"Variance"]), rho=1/r, bet=fixef(fit))
      h.hat <- U %*% drop(fit$coef$random[[1]])
      Vinv <- chol2inv(chol(In + as.numeric(VarCorr(fit)[1,"Variance"])/(fit$sigma)^2*Kmat0))
      hnew <- drop(as.numeric(VarCorr(fit)[1,"Variance"])/(fit$sigma)^2*Kmat10%*%Vinv%*%(y - X%*%fixef(fit)))
      
      preds[[i]][, j] <- hnew
      resids[[i]][, j] <- stats::resid(fit)
      h.hat.ests[[i]][, j] <- h.hat
      
      if(verbose) message("Completed: variable", i, ", r value ", j)
    }
  }
  
  res <- list(q.seq = q.seq, r.seq = r.seq, Drange = Drange, Znew = Znew.mat, resids = resids, preds = preds, h.hat = h.hat.ests)
}

#' Plot of exposure-response function from univariate KMR fot
#' 
#' Plot the estimated \code{h(z[m])} estimated from frequentist KMR for \code{r[m]} fixed to specific values 
#' 
#' @inheritParams kmbayes
#' @param fits output from \code{\link{InvestigatePrior}}
#' @param which.z which predictors (columns in \code{Z}) to plot
#' @param which.q which q.values to plot; defaults to all possible
#' @param plot.resid whether to plot the data points
#' @param ylim plotting limits for the y-axis
#' @param ... other plotting arguments
#' @export
#' @import graphics
PlotPriorFits <- function(y, X, Z, fits, which.z = NULL, which.q = NULL, plot.resid = TRUE, ylim = NULL, ...) {
  q.seq <- fits$q.seq
  r.seq <- fits$r.seq
  Znew <- fits$Znew
  preds <- fits$preds
  
  if (is.null(which.z)) which.z <- 1:ncol(Z)
  if (is.null(which.q)) which.q <- 1:length(q.seq)
  
  q.seq <- q.seq[which.q]
  r.seq <- r.seq[which.q]
  Znew <- Znew[, which.z]
  preds <- preds[which.z]
  Z <- Z[, which.z]
  
  if (plot.resid) {
    lm0 <- lm(y ~ X)
    res <- resid(lm0) + coef(lm0)["(Intercept)"]
    if (is.null(ylim)) ylim <- range(res)
  }
  
  opar <- par(mfrow=c(1 + length(which.z), length(which.q)), mar=c(4.1, 4.1, 1.6, 1.1))
  on.exit(par(opar), add = TRUE)
  
  for(r in r.seq) {
    fun_plot <- function(x) exp(-r*x^2)
    curve(fun_plot, main=paste0("r = ", format(round(r,2), digits = 2, nsmall = 2)), ylab="correlation", cex.lab=1.5, cex.main=2, ylim=c(0,1), xlim=c(0, max(Z)), xlab=expression(d[ij]), xname = 'x')
  }
  for(i in 1:ncol(Z)) {
    for(j in seq_along(r.seq)) {
      est <- preds[[i]][, j]
      if (is.null(ylim)) ylim <- range(est, na.rm = TRUE)
      plot(0, type = "n", ylim = ylim, ylab = expression(hat(h)), xlab = colnames(Z)[i], cex.lab = 1.5, xlim = range(Znew), ...)
      if (plot.resid) points(Z[, i], res, col="red", pch=19, cex=0.5)
      lines(Znew[, i], est)
    }
  }
}
