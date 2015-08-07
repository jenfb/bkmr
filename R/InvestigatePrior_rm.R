#' @export
#' @import nlme
InvestigatePrior <- function(y, expos, covar, ngrid = 50, q.seq = c(2, 1, 1/2, 1/4, 1/8, 1/16), verbose = FALSE) {

    zranges <- diff(apply(expos, 2, range))
    Drange <- max(zranges)
    # q.seq <- c(1/2, 1/5, 1/10, 1/20)
    # q.seq <- c(2, 1, 1/2, 1/4)
    # q.seq <- c(2, 1, 1/2, 1/4, 1/6, 1/8, 1/10)
    # q.seq <- c(2, 1, 1/2, 1/4, 1/8, 1/16)
    r.seq <- -log(1 - 0.50)/(q.seq * Drange)^2

    Znew.mat <- matrix(NA, ngrid, ncol(expos))
    preds <- vector("list", ncol(expos))
    resids <- vector("list", ncol(expos))
    h.hat.ests <- vector("list", ncol(expos))
    for(i in 1:ncol(expos)) {
        preds[[i]] <- matrix(NA, ngrid, length(q.seq))
        resids[[i]] <- matrix(NA, nrow(expos), length(q.seq))
        h.hat.ests[[i]] <- matrix(NA, nrow(expos), length(q.seq))
    }
    for(i in 1:ncol(expos)) {
        Zi <- cbind(expos[, i])
        Znew <- cbind(seq(min(Zi), max(Zi), length.out = ngrid))
        Znew.mat[, i] <- Znew
        n0 <- nrow(Zi)
        In <- diag(1,n0,n0)
        n1 <- nrow(Znew)
        Zall <- rbind(Zi, Znew)
        nall <- n0+n1
        for(j in seq_along(r.seq)) {
            r <- r.seq[j]
            Kpart <- as.matrix(dist(Zall))^2
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
            fit <- lme(y ~ -1+covar, random = list(group = pdIdent(~ -1+U)))
            #data.frame(sig=(fit$sigma)^2, tau=as.numeric(VarCorr(fit)[1,"Variance"]), rho=1/r, bet=fixef(fit))
            h.hat <- U %*% drop(fit$coef$random[[1]])
            Vinv <- chol2inv(chol(In + as.numeric(VarCorr(fit)[1,"Variance"])/(fit$sigma)^2*Kmat0))
            hnew <- drop(as.numeric(VarCorr(fit)[1,"Variance"])/(fit$sigma)^2*Kmat10%*%Vinv%*%(y - covar%*%fixef(fit)))

            preds[[i]][, j] <- hnew
            resids[[i]][, j] <- resid(fit)
            h.hat.ests[[i]][, j] <- h.hat

            if(verbose) message("Completed: variable", i, ", r value ", j)
        }
    }

    res <- list(q.seq = q.seq, r.seq = r.seq, Drange = Drange, Znew = Znew.mat, resids = resids, preds = preds, h.hat = h.hat.ests)
}

#' PlotPriorFits
#'
#' @param y
#' @param covar
#' @param expos
#' @param fits
#' @param which.expos
#' @param which.q
#' @param plot.resid
#' @param ylim
#' @param ...
#'
#' @return fitted objects
#' @export
#'
PlotPriorFits <- function(y, covar, expos, fits, which.expos = NULL, which.q = NULL, plot.resid = TRUE, ylim = NULL, ...) {
    q.seq <- fits$q.seq
    r.seq <- fits$r.seq
    Znew <- fits$Znew
    preds <- fits$preds

    if (is.null(which.expos)) which.expos <- 1:ncol(expos)
    if (is.null(which.q)) which.q <- 1:length(q.seq)

    q.seq <- q.seq[which.q]
    r.seq <- r.seq[which.q]
    Znew <- Znew[, which.expos]
    preds <- preds[which.expos]
    expos <- expos[, which.expos]

    if (plot.resid) {
        lm0 <- lm(y ~ covar)
        res <- resid(lm0) + coef(lm0)["(Intercept)"]
        if (is.null(ylim)) ylim <- range(res)
    }

    opar <- par(mfrow=c(1 + length(which.expos), length(which.q)), mar=c(4.1, 4.1, 1.6, 1.1))
    on.exit(par(opar), add = TRUE)

    for(r in r.seq) {
        curve(exp(-r*x^2), main=round(r,2), ylab="correlation", cex.lab=1.5, cex.main=2, ylim=c(0,1), xlim=c(0, max(expos)), xlab=expression(d[ij]))
    }
    for(i in 1:ncol(expos)) {
        for(j in seq_along(r.seq)) {
            est <- preds[[i]][, j]
            if (is.null(ylim)) ylim <- range(est, na.rm = TRUE)
            plot(0, type = "n", ylim = ylim, ylab = expression(hat(h)), xlab = colnames(expos)[i], cex.lab = 1.5, xlim = range(Znew), ...)
            if (plot.resid) points(expos[, i], res, col="red", pch=19, cex=0.5)
            lines(Znew[, i], est)
        }
    }
}
