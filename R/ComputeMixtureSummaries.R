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
OverallRiskSummaries <- function(fit, y, expos, covar, qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, preds.method = "approx") {
    point1 <- apply(expos, 2, quantile, q.fixed)
    if(preds.method == "approx") {
        preds.fun <- function(znew) ComputePostmeanHnew(fit = fit, y = y, expos = expos, covar = covar, exposNew = znew)
        riskSummary <- riskSummary.approx
    } else if(preds.method == "samp") {
        preds.fun <- function(znew, ...) SampleHnew(Znew = znew, fit = fit, expos = expos, covar = covar, y = y, ...)
        riskSummary <- riskSummary.samp
    }
    risks.overall <- t(sapply(qs, function(quant) riskSummary(point1 = point1, point2 = apply(expos, 2, quantile, quant), preds.fun = preds.fun)))
    risks.overall <- data.frame(quantile = qs, risks.overall)
}

SingPolRiskSummary <- function(whichpol = 1, fit, expos, covar, y, qs.diff = c(0.25, 0.75), q.fixed = 0.5, preds.method = "approx", ...) {
    point2 <- point1 <- apply(expos, 2, quantile, q.fixed)
    point2[whichpol] <- quantile(expos[, whichpol], qs.diff[2])
    point1[whichpol] <- quantile(expos[, whichpol], qs.diff[1])
    # point1 <- makePoint(whichpol, expos, qs.diff[1], q.fixed)
    # point2 <- makePoint(whichpol, expos, qs.diff[2], q.fixed)
    if(preds.method == "approx") {
        preds.fun <- function(znew) ComputePostmeanHnew(fit = fit, X = covar, Z = expos, y = y, Znew = znew)
        riskSummary <- riskSummary.approx
    } else if(preds.method == "samp") {
        preds.fun <- function(znew, ...) SampleHnew(Znew = znew, fit = fit, expos = expos, covar = covar, y = y, ...)
        riskSummary <- riskSummary.samp
    }
    riskSummary(point1 = point1, point2 = point2, preds.fun = preds.fun, ...)
}

SingPolRiskSummaries <- function(fit, expos, covar, y, pollutants = 1:ncol(expos), qs.diff = c(0.25, 0.75), q.fixed = c(0.25, 0.50, 0.75), preds.method = "approx", expos.names = colnames(expos), ...) {
    require(dplyr)

    if(is.null(expos.names)) expos.names <- paste0("expos", 1:ncol(expos))

    df <- data_frame()
    for(i in seq_along(q.fixed)) {
        for(j in seq_along(pollutants)) {
            risk <- SingPolRiskSummary(whichpol = pollutants[j], fit = fit, expos = expos, covar = covar, y = y, qs.diff = qs.diff, q.fixed = q.fixed[i], preds.method = preds.method, ...)
            df0 <- data_frame(q.fixed = q.fixed[i], exposure = expos.names[j], est = risk["est"], se = risk["se"])
            df <- bind_rows(df, df0)
        }
    }
    df <- mutate(df, exposure = factor(exposure, levels = expos.names[pollutants]), q.fixed = as.factor(q.fixed))
    attr(df, "qs.diff") <- qs.diff
    df
}





SingPolIntSummary <- function(whichpol = 1, fit, expos, covar, y, qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75), preds.method = "approx", ...) {
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
        preds.fun <- function(znew) ComputePostmeanHnew(fit = fit, X = covar, Z = expos, y = y, Znew = znew)
        interactionSummary <- interactionSummary.approx
    } else if(preds.method == "samp") {
        preds.fun <- function(znew, ...) SampleHnew(Znew = znew, fit = fit, expos = expos, covar = covar, y = y, ...)
        interactionSummary <- interactionSummary.samp
    }
    interactionSummary(newz.q1, newz.q2, preds.fun, ...)
}

SingPolIntSummaries <- function(fit, expos, covar, y, pollutants = 1:ncol(expos), qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75), preds.method = "approx", expos.names = colnames(expos), ...) {
    require(dplyr)

    if(is.null(expos.names)) expos.names <- paste0("expos", 1:ncol(expos))

    ints <- sapply(pollutants, function(whichpol)
        SingPolIntSummary(whichpol = whichpol, fit = fit, expos = expos, covar = covar, y = y, qs.diff = qs.diff, qs.fixed = qs.fixed, preds.method, ...)
    )

    df <- data_frame(exposure = factor(expos.names[pollutants], levels = expos.names), est = ints["est", ], se = ints["se", ])
}






