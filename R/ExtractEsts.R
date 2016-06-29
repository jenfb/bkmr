#' Compute summary statistics
#'
#' @param q vector of quantiles
#' @param s vector of posterior samples
SummarySamps <- function(s, q = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
    qs <- quantile(s, q)
    names(qs) <- paste0("q_", 100*q)
    summ <- c(mean = mean(s), sd = sd(s), qs)
    summ <- matrix(summ, nrow = 1, dimnames = list(NULL, names(summ)))
}

#' Extract summary statistics
#'
#' Obtain summary statistics of each parameter from the BKMR fit
#'
#' @param fit An object containing the results returned by a the \code{kmbayes} function 
#' @param q vector of quantiles
#' @param sel logical expression indicating samples to keep; defaults to keeping the second half of all samples 
#'
#' @export
ExtractEsts <- function(fit, q = c(0.025, 0.25, 0.5, 0.75, 0.975), sel = NULL) {
  if (inherits(fit, "bkmrfit")) {
    if (is.null(sel)) {
      sel <- with(fit, seq(floor(iter/2) + 1, iter))
    }
    sigsq.eps <- SummarySamps(fit$sigsq.eps[sel], q = q)
    rownames(sigsq.eps) <- "sigsq.eps"
    
    r <- t(apply(fit$r[sel, ], 2, SummarySamps, q = q))
    rownames(r) <- paste0("r", 1:nrow(r))
    
    beta <- t(apply(fit$beta[sel, , drop = FALSE], 2, SummarySamps, q = q))
    
    lambda <- t(apply(fit$lambda[sel, ,drop = FALSE], 2, SummarySamps, q = q))
    if (nrow(lambda) > 1) {
      rownames(lambda) <- paste0("lambda", 1:nrow(lambda))
    } else {
      rownames(lambda) <- "lambda"
    }
    
    h <- t(apply(fit$h.hat[sel, ], 2, SummarySamps, q = q))
    rownames(h) <- paste0("h", 1:nrow(h))
    
    if (!is.null(fit$hnew)) {
      hnew <- t(apply(fit$hnew[sel, ], 2, SummarySamps, q = q))
      rownames(hnew) <- paste0("hnew", 1:nrow(hnew))
    }
  }
  
  if (nrow(beta) > 1) {
    rownames(beta) <- paste0("beta", 1:nrow(beta))
  } else {
    rownames(beta) <- "beta"
  }
  
  colnames(beta) <- colnames(sigsq.eps)
  colnames(r) <- colnames(sigsq.eps)
  colnames(lambda) <- colnames(sigsq.eps)
  colnames(h) <- colnames(sigsq.eps)
  if (!is.null(fit$hnew)) {
    colnames(hnew) <- colnames(sigsq.eps)
  }
  
  ret <- list(sigsq.eps = data.frame(sigsq.eps), beta = beta, lambda = lambda, h = h, r = r)
  if (!is.null(fit$hnew)) ret$hnew <- hnew
  ret
}

#' Extract samples
#'
#' Extract samples of each parameter from the BKMR fit
#'
#' @inheritParams ExtractEsts
#'
#' @export
ExtractSamps <- function(fit, sel = NULL) {
  if (inherits(fit, "bkmrfit")) {
    if (is.null(sel)) {
      sel <- with(fit, seq(floor(iter/2) + 1, iter))
    }
    
    sigsq.eps <- fit$sigsq.eps[sel]
    sig.eps <- sqrt(sigsq.eps)
    r <- fit$r[sel, ]
    beta <- fit$beta[sel, ]
    lambda <- fit$lambda[sel, ]
    tau <- lambda*sigsq.eps
    h <- fit$h.hat[sel, ]
    if (!is.null(fit$hnew)) hnew <- fit$hnew[sel, ]
  }
  
    if (!is.null(ncol(beta))) colnames(beta) <- paste0("beta", 1:ncol(beta))
    colnames(r) <- paste0("r", 1:ncol(r))
    colnames(h) <- paste0("h", 1:ncol(h))
    if (!is.null(fit$hnew)) colnames(hnew) <- paste0("hnew", 1:ncol(hnew))

    res <- list(sigsq.eps = sigsq.eps, sig.eps = sig.eps, r = r, beta = beta, lambda = lambda, tau = tau, h = h)
    if (!is.null(fit$hnew)) res$hnew <- hnew
    res
}
