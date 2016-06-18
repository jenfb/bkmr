#' Calculate group-specific posterior inclusion probabilities
#'
#' Calculate posterior inclusion probabilities for each group of variables
#' 
#' @inheritParams ExtractEsts
#' @export
CalcGroupPIPs <- function(fit, sel = NULL) {
    groups <- fit$groups
    if (is.null(groups)) {
        stop("Cannot compute group-specific posterior inclusion probabilities; BKMR was not run with variable groups")
    }
    if (is.null(sel)) {
        sel <- with(fit, seq(floor(iter/2) + 1, iter))
    }
    grps <- unique(groups)
    groupincl.probs <- sapply(grps, function(x) mean(rowSums(fit$delta[sel, groups == x, drop = FALSE]) > 0))
    groupincl.probs
}

#' Calculate conditional predictor specific posterior inclusion probabilities
#'
#' For those predictors within a multi-preditor group, as defined using the \code{groups} argument, the posterior inclusion probabilities for the predictor conditional on the group being selected into the model.
#' 
#' @inheritParams ExtractEsts
#' @export
CalcWithinGroupPIPs <- function(fit, sel = NULL) {
    groups <- fit$groups
    if (is.null(groups)) {
        stop("Cannot compute group-specific posterior inclusion probabilities; BKMR was not run with variable groups")
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


#' Calculate variable-specific posterior inclusion probabilities
#' 
#' Calculate variable-specific posterior inclusion probabilities from BKMR model fit
#' 
#' @inheritParams ExtractEsts
#' @export
CalcPIPs <- function(fit, sel = NULL) {
  if (inherits(fit, "bkmrfit")) {
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
#' @inheritParams ExtractEsts
#' @param z.names optional argument providing the names of the variables included in the \code{h} function.
#'
#' @return a data frame with the variable-specific PIPs for BKMR fit with component-wise variable selection, and with the group-specific and conditional (within-group) PIPs for BKMR fit with hierarchical variable selection.
#' @export
ExtractPIPs <- function(fit, sel = NULL, z.names = NULL) {
  if (inherits(fit, "bkmrfit")) {
    if (!fit$varsel) {
      stop("This model was not fit with variable selection.")
    }
    if (is.null(sel)) {
      sel <- with(fit, seq(floor(iter/2) + 1, iter))
    }
    if (is.null(z.names)) {
      z.names <- colnames(fit$Z)
    }
    if (is.null(z.names)) {
      z.names <- paste0("z", 1:ncol(fit$Z))
    }
    df <- data.frame(variable = z.names, stringsAsFactors = FALSE)
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

