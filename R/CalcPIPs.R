#' Calculate group-specific posterior inclusion probabilities
#'
#' Calculate posterior inclusion probabilities for each group of pollutants
#'
#' @param fit
#' @param sel
#'
#' @export
CalcGroupPIPs <- function(fit, sel = NULL) {
    groups <- fit$groups
    if (is.null(groups)) {
        stop("Cannot compute group-specific posterior inclusion probabilities; BKMR was not run with pollutant groups")
    }
    if (is.null(sel)) {
        sel <- with(fit, seq(floor(iter/2) + 1, iter))
    }
    grps <- unique(groups)
    groupincl.probs <- sapply(grps, function(x) mean(rowSums(fit$delta[sel, groups == x, drop = FALSE]) > 0))
    groupincl.probs
}

CalcWithinGroupPIPs <- function(fit, sel = NULL) {
    groups <- fit$groups
    if (is.null(groups)) {
        stop("Cannot compute group-specific posterior inclusion probabilities; BKMR was not run with pollutant groups")
    }
    if (is.null(sel)) {
        sel <- with(fit, seq(floor(iter/2) + 1, iter))
    }

    # grps <- unique(groups)
    # ngrps <- sapply(grps, function(x) sum(groups == x))
    condprobs.group <- rep(NA, length(groups))
    for (i in unique(groups)) {
        condprobs.group[groups == i] <- colMeans(fit$delta[sel, ][rowSums(fit$delta[sel, groups == i, drop = FALSE]) > 0, groups == i, drop = FALSE])
    }

    condprobs.group
}


#' Calculate pollutant-specific posterior inclusion probabilities
#'
#' @param fit
#' @param sel
#'
#' @export
#'
CalcPIPs <- function(fit, sel = NULL) {
    if (inherits(fit, "stanfit")) {
        stop("This model was fit using rstan. Stan does not currently support discrete parameter space. Therefore posterior inclusion probabilities are not directly computed.")
    } else if (inherits(fit, "bkmrfit")) {
        if (is.null(sel)) {
            sel <- with(fit, seq(floor(iter/2) + 1, iter))
        }
        groups <- fit$groups
        if (is.null(groups)) {
            ret <- colMeans(fit$delta[sel, , drop = FALSE])
        }

    }
    ret
}

#' Extract posterior inclusion probabilities (PIPs) from BKMR model fit
#'
#' Extract posterior inclusion probabilities (PIPs) from Bayesian Kernel Machine Regression (BKMR) model fit
#'
#' @param fit fitted BKMR model with either component-wise or hierarchical variable selection
#' @param sel optional argument selecting which iterations of the MCMC sampler to keep
#' @param expos.names optional argument providing the names of the exposure variables included in the \code{h} function.
#'
#' @return a data frame with the pollutant-specific PIPs for BKMR fit with component-wise variable selection, and with the group-specific and conditoinal (within-group) PIPs for BKMR fit with hierarchical variable selection.
#' @export
ExtractPIPs <- function(fit, sel = NULL, expos.names = NULL) {
    if (inherits(fit, "stanfit")) {
        stop("This model was fit using rstan. Stan does not currently support discrete parameter space. Therefore posterior inclusion probabilities are not directly computed.")
    } else if (inherits(fit, "bkmrfit")) {
        if (!fit$modsel) {
            stop("This model was not fit with variable selection.")
        }
        if (is.null(sel)) {
            sel <- with(fit, seq(floor(iter/2) + 1, iter))
        }
        if (is.null(expos.names)) {
            expos.names <- colnames(fit$Z)
        }
        if (is.null(expos.names)) {
            expos.names <- paste0("expos", 1:ncol(expos))
        }
        df <- data.frame(exposure = expos.names, stringsAsFactors = FALSE)
        groups <- fit$groups
        if (is.null(groups)) {
            df$PIP <- colMeans(fit$delta[sel, , drop = FALSE])
        } else {

            ## group-specific posterior inclusion probability
            df$group <- groups
            grps <- unique(groups)
            groupincl.probs <- sapply(grps, function(x) mean(rowSums(fit$delta[sel, groups == x, drop = FALSE]) > 0))
            df.group <- dplyr::data_frame(group = grps,
                                          groupPIP = groupincl.probs)
            df <- dplyr::inner_join(df, df.group, by = "group")

            ## within-group conditional PIP
            condprobs.group <- rep(NA, length(groups))
            for (i in unique(groups)) {
                condprobs.group[groups == i] <- colMeans(fit$delta[sel, ][rowSums(fit$delta[sel, groups == i, drop = FALSE]) > 0, groups == i, drop = FALSE])
            }
            df$condPIP <- condprobs.group
        }
    }
    df
}

