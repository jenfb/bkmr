#' @param s vector of posterior samples
SummarySamps <- function(s, q = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
    qs <- quantile(s, q)
    names(qs) <- paste0("q_", 100*q)
    summ <- c(mean = mean(s), sd = sd(s), qs)
    summ <- matrix(summ, nrow = 1, dimnames = list(NULL, names(summ)))
}

#' Obtain summary statistics of each parameter from the BKMR fit
#'
#' @export
ExtractEsts <- function(fit, q = c(0.025, 0.25, 0.5, 0.75, 0.975), sel = NULL) {
    if(inherits(fit, "stanfit")) {
        samps <- rstan::extract(fit)

        sigsq.eps <- SummarySamps(samps$sigma_sq, q = q)
        rownames(sigsq.eps) <- "sigsq.eps"

        r <- t(apply(samps$r, 2, SummarySamps, q = q))
        rownames(r) <- paste0("r", 1:nrow(r))

        beta <- t(apply(samps$beta, 2, SummarySamps, q = q))

        lambda <- SummarySamps(samps$tau/samps$sigma_sq, q = q)
        rownames(lambda) <- "lambda"

        h <- t(apply(samps$h, 2, SummarySamps, q = q))
        rownames(h) <- paste0("h", 1:nrow(h))

    } else if (inherits(fit, "bkmrfit")) {
        if (is.null(sel)) {
            sel <- with(fit, seq(floor(iter/2) + 1, iter))
        }
        sigsq.eps <- SummarySamps(fit$sigsq.eps[sel], q = q)
        rownames(sigsq.eps) <- "sigsq.eps"

        r <- t(apply(fit$r[sel, ], 2, SummarySamps, q = q))
        rownames(r) <- paste0("r", 1:nrow(r))

        beta <- t(apply(fit$beta[sel, , drop=FALSE], 2, SummarySamps, q = q))

        lambda <- t(apply(fit$lambda[sel, ,drop = FALSE], 2, SummarySamps, q = q))
        if (nrow(lambda) > 1) {
            rownames(lambda) <- paste0("lambda", 1:nrow(lambda))
        } else {
            rownames(lambda) <- "lambda"
        }

        h <- t(apply(fit$h[sel, ], 2, SummarySamps, q = q))
        rownames(h) <- paste0("h", 1:nrow(h))
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

    list(sigsq.eps = data.frame(sigsq.eps), beta = beta, lambda = lambda, h = h, r = r)
}

#' Extract samples of each parameter from the BKMR fit
#'
#' @export
ExtractSamps <- function(fit, sel = NULL) {
    if(inherits(fit, "stanfit")) {
        samps <- rstan::extract(fit)

        sig.eps <- samps$sigma
        sigsq.eps <- samps$sigma_sq
        r <- samps$r
        beta <- samps$beta
        tau <- samps$tau
        lambda <- samps$tau/samps$sigma_sq
        h <- samps$h

    } else if (inherits(fit, "bkmrfit")) {
        if (is.null(sel)) {
            sel <- with(fit, seq(floor(iter/2) + 1, iter))
        }

        sigsq.eps <- fit$sigsq.eps[sel]
        sig.eps <- sqrt(sigsq.eps)
        r <- fit$r[sel, ]
        beta <- fit$beta[sel, ]
        lambda <- fit$lambda[sel, ]
        tau <- lambda*sigsq.eps
        h <- fit$h[sel, ]
    }

    if(!is.null(ncol(beta))) colnames(beta) <- paste0("beta", 1:ncol(beta))
    colnames(r) <- paste0("r", 1:ncol(r))
    colnames(h) <- paste0("h", 1:ncol(h))

    res <- list(sigsq.eps = sigsq.eps, sig.eps = sig.eps, r = r, beta = beta, lambda = lambda, tau = tau, h = h)
}
