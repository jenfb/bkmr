riskSummary.approx <- function(point1, point2, preds.fun, ...) {
    cc <- c(-1, 1)
    newz <- rbind(point1, point2)
    preds <- preds.fun(newz, ...)
    diff <- drop(cc %*% preds$postmean)
    diff.sd <- drop(sqrt(cc %*% preds$postvar %*% cc))
    c(est = diff, sd = diff.sd)
}
riskSummary.samp <- function(point1, point2, preds.fun, ...) {
    cc <- c(-1, 1)
    newz <- rbind(point1, point2)
    preds <- preds.fun(newz, ...)
    diff.preds <- drop(preds %*% cc)
    c(est = mean(diff.preds), sd = sd(diff.preds))
}
interactionSummary.approx <- function(newz.q1, newz.q2, preds.fun, ...) {
    cc <- c(-1*c(-1, 1), c(-1, 1))
    newz <- rbind(newz.q1, newz.q2)
    preds <- preds.fun(newz, ...)
    int <- drop(cc %*% preds$postmean)
    int.se <- drop(sqrt(cc %*% preds$postvar %*% cc))
    c(est = int, sd = int.se)
}
interactionSummary.samp <- function(newz.q1, newz.q2, preds.fun, ...) {
    cc <- c(-1*c(-1, 1), c(-1, 1))
    newz <- rbind(newz.q1, newz.q2)
    preds <- preds.fun(newz, ...)
    int.preds <- drop(preds %*% cc)
    c(est = mean(int.preds), sd = sd(int.preds))
}



#' OverallRiskSummaries
#' 
#' Compare estimated \code{h} function when all predictors are at a particular percentile to when all are at the 50th percentile
#' @inheritParams kmbayes
#' @inheritParams ExtractEsts
#'
#' @param fit 
#' @param y 
#' @param Z 
#' @param X 
#' @param qs 
#' @param q.fixed 
#' @param preds.method 
#' @param sel 
#'
#' @export
OverallRiskSummaries <- function(fit, y, Z, X, qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, preds.method = "approx", sel = NULL) {
    point1 <- apply(Z, 2, quantile, q.fixed)
    if(preds.method == "approx") {
        preds.fun <- function(znew) ComputePostmeanHnew(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel)
        riskSummary <- riskSummary.approx
    } else if(preds.method == "samp") {
        stop("not yet implemented")
#         preds.fun <- function(znew, ...) SampleHnew(Znew = znew, fit = fit, Z = Z, X = X, y = y, ...)
#         riskSummary <- riskSummary.samp
    }
    risks.overall <- t(sapply(qs, function(quant) riskSummary(point1 = point1, point2 = apply(Z, 2, quantile, quant), preds.fun = preds.fun)))
    risks.overall <- data.frame(quantile = qs, risks.overall)
}

#' VarRiskSummary
#' 
#' Compare estimated \code{h} function when a single variable (or a set of variables) is at the 75th versus 25th percentile, when all of the other variables are fixed at a particular percentile
#' @inheritParams kmbayes
#'
#' @export
#'
#' @param fit 
#' @param y 
#' @param Z 
#' @param X 
#' @param qs.diff 
#' @param q.fixed 
#' @param preds.method 
#' @param sel 
#' @param ... 
#' @param whichz a scalar or vector selecting which Z variables to compute the summary for (the other variables in Z will be fixed at the value \code{q.fixed})
VarRiskSummary <- function(whichz = 1, fit, y, Z, X, qs.diff = c(0.25, 0.75), q.fixed = 0.5, preds.method = "approx", sel = NULL, ...) {
    point2 <- point1 <- apply(Z, 2, quantile, q.fixed)
    point2[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[2])
    point1[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[1])
    # point1 <- makePoint(whichz, Z, qs.diff[1], q.fixed)
    # point2 <- makePoint(whichz, Z, qs.diff[2], q.fixed)
    if(preds.method == "approx") {
        preds.fun <- function(znew) ComputePostmeanHnew(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel)
        riskSummary <- riskSummary.approx
    } else if(preds.method == "samp") {
        stop("not yet implemented")
#         preds.fun <- function(znew, ...) SampleHnew(Znew = znew, fit = fit, Z = Z, X = X, y = y, ...)
#         riskSummary <- riskSummary.samp
    }
    riskSummary(point1 = point1, point2 = point2, preds.fun = preds.fun, ...)
}

#' SingVarRiskSummaries
#' 
#' @inheritParams kmbayes
#' @param fit 
#'
#' @param y 
#' @param Z 
#' @param X 
#' @param which.z 
#' @param qs.diff 
#' @param q.fixed 
#' @param preds.method 
#' @param sel 
#' @param z.names 
#' @param ... 
#'
#' @export
SingVarRiskSummaries <- function(fit, y, Z, X, which.z = 1:ncol(Z), qs.diff = c(0.25, 0.75), q.fixed = c(0.25, 0.50, 0.75), preds.method = "approx", sel = NULL, z.names = colnames(Z), ...) {
    if(is.null(z.names)) z.names <- paste0("z", 1:ncol(Z))

    df <- dplyr::data_frame()
    for(i in seq_along(q.fixed)) {
        for(j in seq_along(which.z)) {
            risk <- VarRiskSummary(whichz = which.z[j], fit = fit, y = y, Z = Z, X = X, qs.diff = qs.diff, q.fixed = q.fixed[i], preds.method = preds.method, sel = sel, ...)
            df0 <- dplyr::data_frame(q.fixed = q.fixed[i], variable = z.names[j], est = risk["est"], sd = risk["sd"])
            df <- dplyr::bind_rows(df, df0)
        }
    }
    df <- dplyr::mutate_(df, variable = ~factor(variable, levels = z.names[which.z]), q.fixed = ~as.factor(q.fixed))
    attr(df, "qs.diff") <- qs.diff
    df
}




#' SingVarIntSummary
#' 
#' Compare the single-predictor health risks when all of the other predictors in Z are fixed to their 75th percentile to when all of the other predictors in Z are fixed to their 25th percentile.
#' @inheritParams kmbayes
#'
#' @param whichz 
#' @param fit 
#' @param y 
#' @param Z 
#' @param X 
#' @param qs.diff 
#' @param qs.fixed 
#' @param preds.method 
#' @param sel 
#' @param ... 
#'
#' @export
SingVarIntSummary <- function(whichz = 1, fit, y, Z, X, qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75), preds.method = "approx", sel = NULL, ...) {
    q.fixed <- qs.fixed[1]
    point2 <- point1 <- apply(Z, 2, quantile, q.fixed)
    point2[whichz] <- quantile(Z[, whichz], qs.diff[2])
    point1[whichz] <- quantile(Z[, whichz], qs.diff[1])
    newz.q1 <- rbind(point1, point2)

    q.fixed <- qs.fixed[2]
    point2 <- point1 <- apply(Z, 2, quantile, q.fixed)
    point2[whichz] <- quantile(Z[, whichz], qs.diff[2])
    point1[whichz] <- quantile(Z[, whichz], qs.diff[1])
    newz.q2 <- rbind(point1, point2)

    if(preds.method == "approx") {
        preds.fun <- function(znew) ComputePostmeanHnew(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel)
        interactionSummary <- interactionSummary.approx
    } else if(preds.method == "samp") {
        stop("not yet implemented")
#         preds.fun <- function(znew, ...) SampleHnew(Znew = znew, fit = fit, Z = Z, X = X, y = y, ...)
#         interactionSummary <- interactionSummary.samp
    }
    interactionSummary(newz.q1, newz.q2, preds.fun, ...)
}

#' SingVarIntSummaries
#' 
#' @inheritParams kmbayes
#' @param fit 
#'
#' @param y 
#' @param Z 
#' @param X 
#' @param which.z 
#' @param qs.diff 
#' @param qs.fixed 
#' @param preds.method 
#' @param sel 
#' @param z.names 
#' @param ... 
#'
#' @export
SingVarIntSummaries <- function(fit, y, Z, X, which.z = 1:ncol(Z), qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75), preds.method = "approx", sel = NULL, z.names = colnames(Z), ...) {
    if(is.null(z.names)) z.names <- paste0("z", 1:ncol(Z))

    ints <- sapply(which.z, function(whichz)
        SingVarIntSummary(whichz = whichz, fit = fit, Z = Z, X = X, y = y, qs.diff = qs.diff, qs.fixed = qs.fixed, preds.method, sel = sel, ...)
    )

    df <- dplyr::data_frame(variable = factor(z.names[which.z], levels = z.names), est = ints["est", ], sd = ints["sd", ])
}



