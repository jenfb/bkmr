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



#' Calculate overall risk summaries
#' 
#' Compare estimated \code{h} function when all predictors are at a particular quantile to when all are at a second fixed quantile
#' @inheritParams kmbayes
#' @inheritParams ComputePostmeanHnew
#' @inherit ComputePostmeanHnew details
#' @param qs vector of quantiles at which to calculate the overall risk summary 
#' @param q.fixed a second quantile at which to compare the estimated \code{h} function
#' @export
#' @return a data frame containing the (posterior mean) estimate and posterior standard deviation of the overall risk measures
#' @examples
#' ## First generate dataset
#' set.seed(111)
#' dat <- SimData(n = 50, M = 4)
#' y <- dat$y
#' Z <- dat$Z
#' X <- dat$X
#' 
#' ## Fit model with component-wise variable selection
#' ## Using only 100 iterations to make example run quickly
#' ## Typically should use a large number of iterations for inference
#' set.seed(111)
#' fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE)
#' 
#' risks.overall <- OverallRiskSummaries(fit = fitkm, qs = seq(0.25, 0.75, by = 0.05), 
#' q.fixed = 0.5, method = "exact")
OverallRiskSummaries <- function(fit, y = NULL, Z = NULL, X = NULL, qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, method = "approx", sel = NULL) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
  }
  
  point1 <- apply(Z, 2, quantile, q.fixed)
  if (method %in% c("approx", "exact")) {
    preds.fun <- function(znew) ComputePostmeanHnew(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel, method = method)
    riskSummary <- riskSummary.approx
  } else {
    stop("method must be one of c('approx', 'exact')")
  }
  risks.overall <- t(sapply(qs, function(quant) riskSummary(point1 = point1, point2 = apply(Z, 2, quantile, quant), preds.fun = preds.fun)))
  risks.overall <- data.frame(quantile = qs, risks.overall)
}

#Compare estimated \code{h} function when a single variable (or a set of variables) is at the 75th versus 25th percentile, when all of the other variables are fixed at a particular percentile
VarRiskSummary <- function(whichz = 1, fit, y = NULL, Z = NULL, X = NULL, qs.diff = c(0.25, 0.75), q.fixed = 0.5, method = "approx", sel = NULL, ...) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
  }
  
  point2 <- point1 <- apply(Z, 2, quantile, q.fixed)
  point2[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[2])
  point1[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[1])
  # point1 <- makePoint(whichz, Z, qs.diff[1], q.fixed)
  # point2 <- makePoint(whichz, Z, qs.diff[2], q.fixed)
  if (method %in% c("approx", "exact")) {
    preds.fun <- function(znew) ComputePostmeanHnew(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel, method = method)
    riskSummary <- riskSummary.approx
  }  else {
    stop("method must be one of c('approx', 'exact')")
  }
  riskSummary(point1 = point1, point2 = point2, preds.fun = preds.fun, ...)
}

#' Single Variable Risk Summaries
#' 
#' Compute summaries of the risks associated with a change in a single variable in \code{Z} from a single level (quantile) to a second level (quantile), for the other variables in \code{Z} fixed to a specific level (quantile)
#' 
#' @inheritParams kmbayes
#' @inheritParams ExtractEsts
#' @inheritParams OverallRiskSummaries
#' @inherit ComputePostmeanHnew details
#' @param qs.diff vector indicating the two quantiles \code{q_1} and \code{q_2} at which to compute \code{h(z_{q2}) - h(z_{q1})}
#' @param q.fixed vector of quantiles at which to fix the remaining predictors in \code{Z}
#' @param z.names optional vector of names for the columns of \code{z}
#' @param ... other arguments to pass on to the prediction function
#' @param which.z vector indicating which variables (columns of \code{Z}) for which the summary should be computed
#' @export
#' 
#' @return a data frame containing the (posterior mean) estimate and posterior standard deviation of the single-predictor risk measures
#' 
#' @examples
#' ## First generate dataset
#' set.seed(111)
#' dat <- SimData(n = 50, M = 4)
#' y <- dat$y
#' Z <- dat$Z
#' X <- dat$X
#' 
#' ## Fit model with component-wise variable selection
#' ## Using only 100 iterations to make example run quickly
#' ## Typically should use a large number of iterations for inference
#' set.seed(111)
#' fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE)
#' 
#' risks.singvar <- SingVarRiskSummaries(fit = fitkm, method = "exact")
SingVarRiskSummaries <- function(fit, y = NULL, Z = NULL, X = NULL, which.z = 1:ncol(Z), qs.diff = c(0.25, 0.75), q.fixed = c(0.25, 0.50, 0.75), method = "approx", sel = NULL, z.names = colnames(Z), ...) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
  }
  
  if(is.null(z.names)) z.names <- paste0("z", 1:ncol(Z))
  
  df <- dplyr::tibble()
  for(i in seq_along(q.fixed)) {
    for(j in seq_along(which.z)) {
      risk <- VarRiskSummary(whichz = which.z[j], fit = fit, y = y, Z = Z, X = X, qs.diff = qs.diff, q.fixed = q.fixed[i], method = method, sel = sel, ...)
      df0 <- dplyr::tibble(q.fixed = q.fixed[i], variable = z.names[j], est = risk["est"], sd = risk["sd"])
      df <- dplyr::bind_rows(df, df0)
    }
  }
  #df <- dplyr::mutate_(df, variable = ~factor(variable, levels = z.names[which.z]), q.fixed = ~as.factor(q.fixed))
  df <- dplyr::mutate_at(df, "variable", function(x) factor(x, levels = z.names[which.z]))
  df <- dplyr::mutate_at(df, "q.fixed", function(x) as.factor(x))
  attr(df, "qs.diff") <- qs.diff
  df
}

SingVarIntSummary <- function(whichz = 1, fit, y = NULL, Z = NULL, X = NULL, qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75), method = "approx", sel = NULL, ...) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
  }
  
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
  
  if (method %in% c("approx", "exact")) {
    preds.fun <- function(znew) ComputePostmeanHnew(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel, method = method)
    interactionSummary <- interactionSummary.approx
  } else {
    stop("method must be one of c('approx', 'exact')")
  }
  interactionSummary(newz.q1, newz.q2, preds.fun, ...)
}

#' Single Variable Interaction Summaries
#' 
#' Compare the single-predictor health risks when all of the other predictors in Z are fixed to their a specific quantile to when all of the other predictors in Z are fixed to their a second specific quantile.
#' @inheritParams kmbayes
#' @inheritParams ExtractEsts
#' @inheritParams SingVarRiskSummaries
#' @inherit ComputePostmeanHnew details
#' @param qs.diff vector indicating the two quantiles at which to compute the single-predictor risk summary
#' @param qs.fixed vector indicating the two quantiles at which to fix all of the remaining exposures in \code{Z}
#' @export
#' 
#' @return a data frame containing the (posterior mean) estimate and posterior standard deviation of the single-predictor risk measures
#' 
#' @examples
#' ## First generate dataset
#' set.seed(111)
#' dat <- SimData(n = 50, M = 4)
#' y <- dat$y
#' Z <- dat$Z
#' X <- dat$X
#' 
#' ## Fit model with component-wise variable selection
#' ## Using only 100 iterations to make example run quickly
#' ## Typically should use a large number of iterations for inference
#' set.seed(111)
#' fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE)
#' 
#' risks.int <- SingVarIntSummaries(fit = fitkm, method = "exact")
SingVarIntSummaries <- function(fit, y = NULL, Z = NULL, X = NULL, which.z = 1:ncol(Z), qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75), method = "approx", sel = NULL, z.names = colnames(Z), ...) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
  }
  
  if(is.null(z.names)) z.names <- paste0("z", 1:ncol(Z))
  
  ints <- sapply(which.z, function(whichz)
    SingVarIntSummary(whichz = whichz, fit = fit, Z = Z, X = X, y = y, qs.diff = qs.diff, qs.fixed = qs.fixed, method, sel = sel, ...)
  )
  
  df <- dplyr::tibble(variable = factor(z.names[which.z], levels = z.names), est = ints["est", ], sd = ints["sd", ])
}



