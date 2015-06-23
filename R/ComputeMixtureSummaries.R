riskSummary.approx <- function(point1, point2, preds.fun, ...) {
    cc <- c(-1, 1)
    newz <- rbind(point1, point2)
    preds <- preds.fun(newz, ...)
    diff <- drop(cc %*% preds$postmean)
    diff.se <- drop(sqrt(cc %*% preds$postvar %*% cc))
    c(est = diff, se = diff.se)
}
riskSummary.samp <- function(point1, point2, preds.fun, ...) {
    cc <- c(-1, 1)
    newz <- rbind(point1, point2)
    preds <- preds.fun(newz, ...)
    diff.preds <- drop(preds %*% cc)
    c(est = mean(diff.preds), se = sd(diff.preds))
}
interactionSummary.approx <- function(newz.q1, newz.q2, preds.fun, ...) {
    cc <- c(-1*c(-1, 1), c(-1, 1))
    newz <- rbind(newz.q1, newz.q2)
    preds <- preds.fun(newz, ...)
    int <- drop(cc %*% preds$postmean)
    int.se <- drop(sqrt(cc %*% preds$postvar %*% cc))
    c(est = int, se = int.se)
}
interactionSummary.samp <- function(newz.q1, newz.q2, preds.fun, ...) {
    cc <- c(-1*c(-1, 1), c(-1, 1))
    newz <- rbind(newz.q1, newz.q2)
    preds <- preds.fun(newz, ...)
    int.preds <- drop(preds %*% cc)
    c(est = mean(int.preds), se = sd(int.preds))
}



#' Compare estimated exposure-response function when all exposures are at a particular percentile to when all are at the 50th percentile
#'
#' @export
OverallRiskSummaries <- function(fit, y, expos, covar, qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, preds.method = "approx", sel = NULL) {
    point1 <- apply(expos, 2, quantile, q.fixed)
    if(preds.method == "approx") {
        preds.fun <- function(znew) ComputePostmeanHnew(fit = fit, y = y, expos = expos, covar = covar, exposNew = znew, sel = sel)
        riskSummary <- riskSummary.approx
    } else if(preds.method == "samp") {
        stop("not yet implemented")
#         preds.fun <- function(znew, ...) SampleHnew(Znew = znew, fit = fit, expos = expos, covar = covar, y = y, ...)
#         riskSummary <- riskSummary.samp
    }
    risks.overall <- t(sapply(qs, function(quant) riskSummary(point1 = point1, point2 = apply(expos, 2, quantile, quant), preds.fun = preds.fun)))
    risks.overall <- data.frame(quantile = qs, risks.overall)
}

#' Compare estimated exposure-response function when a single pollutant (or a set of pollutants) is at the 75th versus 25th percentile, when all of the other exposures are fixed at a particular percentile
#'
#' @export
#' @param whichpol a scalar or vector selecting which pollutants to compute the summary for (the other pollutants will be fixed at the value \code{q.fixed})
PolRiskSummary <- function(whichpol = 1, fit, y, expos, covar, qs.diff = c(0.25, 0.75), q.fixed = 0.5, preds.method = "approx", sel = NULL, ...) {
    point2 <- point1 <- apply(expos, 2, quantile, q.fixed)
    point2[whichpol] <- apply(expos[, whichpol, drop = FALSE], 2, quantile, qs.diff[2])
    point1[whichpol] <- apply(expos[, whichpol, drop = FALSE], 2, quantile, qs.diff[1])
    # point1 <- makePoint(whichpol, expos, qs.diff[1], q.fixed)
    # point2 <- makePoint(whichpol, expos, qs.diff[2], q.fixed)
    if(preds.method == "approx") {
        preds.fun <- function(znew) ComputePostmeanHnew(fit = fit, y = y, expos = expos, covar = covar, exposNew = znew, sel = sel)
        riskSummary <- riskSummary.approx
    } else if(preds.method == "samp") {
        stop("not yet implemented")
#         preds.fun <- function(znew, ...) SampleHnew(Znew = znew, fit = fit, expos = expos, covar = covar, y = y, ...)
#         riskSummary <- riskSummary.samp
    }
    riskSummary(point1 = point1, point2 = point2, preds.fun = preds.fun, ...)
}

#' @export
SingPolRiskSummaries <- function(fit, y, expos, covar, pollutants = 1:ncol(expos), qs.diff = c(0.25, 0.75), q.fixed = c(0.25, 0.50, 0.75), preds.method = "approx", sel = NULL, expos.names = colnames(expos), ...) {
    if(is.null(expos.names)) expos.names <- paste0("expos", 1:ncol(expos))

    df <- dplyr::data_frame()
    for(i in seq_along(q.fixed)) {
        for(j in seq_along(pollutants)) {
            risk <- PolRiskSummary(whichpol = pollutants[j], fit = fit, y = y, expos = expos, covar = covar, qs.diff = qs.diff, q.fixed = q.fixed[i], preds.method = preds.method, sel = sel, ...)
            df0 <- dplyr::data_frame(q.fixed = q.fixed[i], exposure = expos.names[j], est = risk["est"], se = risk["se"])
            df <- dplyr::bind_rows(df, df0)
        }
    }
    df <- dplyr::mutate(df, exposure = factor(exposure, levels = expos.names[pollutants]), q.fixed = as.factor(q.fixed))
    attr(df, "qs.diff") <- qs.diff
    df
}




#' Compare the single-pollutant health risks when all of the other pollutants are fixed to their 75th percentile to when all of the other pollutants are fixed to their 25th percentile.
#'
#' @export
SingPolIntSummary <- function(whichpol = 1, fit, y, expos, covar, qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75), preds.method = "approx", sel = NULL, ...) {
    q.fixed <- qs.fixed[1]
    point2 <- point1 <- apply(expos, 2, quantile, q.fixed)
    point2[whichpol] <- quantile(expos[, whichpol], qs.diff[2])
    point1[whichpol] <- quantile(expos[, whichpol], qs.diff[1])
    newz.q1 <- rbind(point1, point2)

    q.fixed <- qs.fixed[2]
    point2 <- point1 <- apply(expos, 2, quantile, q.fixed)
    point2[whichpol] <- quantile(expos[, whichpol], qs.diff[2])
    point1[whichpol] <- quantile(expos[, whichpol], qs.diff[1])
    newz.q2 <- rbind(point1, point2)

    if(preds.method == "approx") {
        preds.fun <- function(znew) ComputePostmeanHnew(fit = fit, y = y, expos = expos, covar = covar, exposNew = znew, sel = sel)
        interactionSummary <- interactionSummary.approx
    } else if(preds.method == "samp") {
        stop("not yet implemented")
#         preds.fun <- function(znew, ...) SampleHnew(Znew = znew, fit = fit, expos = expos, covar = covar, y = y, ...)
#         interactionSummary <- interactionSummary.samp
    }
    interactionSummary(newz.q1, newz.q2, preds.fun, ...)
}

#' @export
SingPolIntSummaries <- function(fit, y, expos, covar, pollutants = 1:ncol(expos), qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75), preds.method = "approx", sel = NULL, expos.names = colnames(expos), ...) {
    if(is.null(expos.names)) expos.names <- paste0("expos", 1:ncol(expos))

    ints <- sapply(pollutants, function(whichpol)
        SingPolIntSummary(whichpol = whichpol, fit = fit, expos = expos, covar = covar, y = y, qs.diff = qs.diff, qs.fixed = qs.fixed, preds.method, sel = sel, ...)
    )

    df <- dplyr::data_frame(exposure = factor(expos.names[pollutants], levels = expos.names), est = ints["est", ], se = ints["se", ])
}


#' Sequentially add pollutants to compute the overall risk, when the set of polutants are all at their 75th vs. 25th percentile, for all of the other pollutants fixed at a particular percentile
#'
#' @export
SeqPolRiskSummaries <- function(fit, y, expos, covar, pollutants = 1:ncol(expos), qs.diff = c(0.25, 0.75), q.fixed = 0.50, preds.method = "approx", sel = NULL, expos.names = colnames(expos), ...) {
    if(is.null(expos.names)) expos.names <- paste0("expos", 1:ncol(expos))

    if (is.character(colnames(expos))) {
        expos.names <- expos.names[match(pollutants, expos.names)]
    } else {
        expos.names <- expos.names[pollutants]
    }

    df <- dplyr::data_frame()
    for(i in seq_along(q.fixed)) {
        for(j in seq_along(pollutants)) {
            risk <- PolRiskSummary(whichpol = pollutants[1:j], fit = fit, y = y, expos = expos, covar = covar, qs.diff = qs.diff, q.fixed = q.fixed[i], preds.method = preds.method, sel = sel, ...)
            df0 <- dplyr::data_frame(q.fixed = q.fixed[i], added.exposure = expos.names[j], est = risk["est"], se = risk["se"])
            df <- dplyr::bind_rows(df, df0)
        }
    }
    df <- dplyr::mutate(df, added.exposure = factor(added.exposure, levels = expos.names), q.fixed = as.factor(q.fixed))
    attr(df, "qs.diff") <- qs.diff
    df
}



