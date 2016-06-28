PredictorResponseUnivarVar <- function(whichz = 1, fit, y, Z, X, preds.method = "approx", ngrid = 50, q.fixed = 0.5, sel = NULL, min.plot.dist = Inf, center = TRUE, z.names = colnames(Z), ...) {

    if (ncol(Z) < 2) stop("requires there to be at least 2 predictor variables")

    if (is.null(z.names)) {
        colnames(Z) <- paste0("z", 1:ncol(Z))
    } else {
        colnames(Z) <- z.names
    }

    ord <- c(whichz, setdiff(1:ncol(Z), whichz))
    z1 <- seq(min(Z[,ord[1]]), max(Z[,ord[1]]), length = ngrid)
    z.others <- lapply(2:ncol(Z), function(x) quantile(Z[,ord[x]], q.fixed))
    z.all <- c(list(z1), z.others)
    newz.grid <- expand.grid(z.all)
    colnames(newz.grid) <- colnames(Z)[ord]
    newz.grid <- newz.grid[,colnames(Z)]

    if (!is.null(min.plot.dist)) {
        mindists <- rep(NA,nrow(newz.grid))
        for (i in seq_along(mindists)) {
            pt <- as.numeric(newz.grid[i, colnames(Z)[ord[1]]])
            dists <- fields::rdist(matrix(pt, nrow = 1), Z[, colnames(Z)[ord[1]]])
            mindists[i] <- min(dists)
        }
    }

    if(preds.method == "approx") {
        preds <- ComputePostmeanHnew(fit = fit, y = y, Z = Z, X = X, Znew = newz.grid, sel = sel)
        preds.plot <- preds$postmean
        se.plot <- sqrt(diag(preds$postvar))
    } else if(preds.method == "samp") {
        stop("not yet implemented")
#         preds <- SampleHnew(Znew = newz.grid, fit = fit, Z = Z, X = X, y = y, ...)
#         preds.plot <- colMeans(preds, na.rm = TRUE)
#         se.plot <- apply(preds, 2, sd, na.rm = TRUE)
    }
    if(center) preds.plot <- preds.plot - mean(preds.plot)
    if(!is.null(min.plot.dist)) {
        preds.plot[mindists > min.plot.dist] <- NA
        se.plot[mindists > min.plot.dist] <- NA
    }

    res <- dplyr::data_frame(z = z1, est = preds.plot, se = se.plot)
}

#' Plot univariate predictor-response function on a new grid of points
#' 
#' Plot univariate predictor-response function on a new grid of points
#' 
#' @inheritParams kmbayes
#' @inheritParams ExtractEsts
#' @inheritParams SingVarRiskSummaries
#' 
#' @param which.z vector identifying which predictors (columns of \code{Z}) should be plotted
#' @param ngrid number of grid points to cover the range of each predictor (column in \code{Z})
#' @param min.plot.dist specifies a minimum distance that a new grid point needs to be from an observed data point in order to compute the prediction; points further than this will not be computed
#' @param center flag for whether to scale the exposure-response function to have mean zero
#'
#' @export
PredictorResponseUnivar <- function(fit, y = NULL, Z = NULL, X = NULL, which.z = 1:ncol(Z), preds.method = "approx", ngrid = 50, q.fixed = 0.5, sel = NULL, min.plot.dist = Inf, center = TRUE, z.names = colnames(Z), ...) {
  
  if (inherits(fit, "bkmrfit")) {
    y <- fit$y
    Z <- fit$Z
    X <- fit$X
  }
 
    df <- dplyr::data_frame()
    for(i in which.z) {
        res <- PredictorResponseUnivarVar(whichz = i, fit = fit, y = y, Z = Z, X = X, preds.method = preds.method, ngrid = ngrid, q.fixed = q.fixed, sel = sel, min.plot.dist = min.plot.dist, center = center, z.names = z.names, ...)
        df0 <- dplyr::mutate(res, variable = z.names[i]) %>%
            dplyr::select_(~variable, ~z, ~est, ~se)
        df <- dplyr::bind_rows(df, df0)
    }
    df$variable <- factor(df$variable, levels = z.names[which.z])
    df
}



PredictorResponseBivarPair <- function(fit, y, Z, X, whichz1 = 1, whichz2 = 2, whichz3, preds.method = "approx", prob = 0.5, q.fixed = 0.5, sel = NULL, ngrid = 50, min.plot.dist = 0.5, center = TRUE, ...) {
    if(ncol(Z) < 3) stop("requires there to be at least 3 Z variables")

    if(is.null(colnames(Z))) colnames(Z) <- paste0("z", 1:ncol(Z))

    if(missing(whichz3)) {
        ord <- c(whichz1, whichz2, setdiff(1:ncol(Z), c(whichz1, whichz2)))
    } else {
        ord <- c(whichz1, whichz2, whichz3, setdiff(1:ncol(Z), c(whichz1, whichz2, whichz3)))
    }
    z1 <- seq(min(Z[,ord[1]]), max(Z[,ord[1]]), length=ngrid)
    z2 <- seq(min(Z[,ord[2]]), max(Z[,ord[2]]), length=ngrid)
    z3 <- quantile(Z[, ord[3]], probs = prob)
    z.all <- c(list(z1), list(z2), list(z3))
    if(ncol(Z) > 3) {
        z.others <- lapply(4:ncol(Z), function(x) quantile(Z[,ord[x]], q.fixed))
        z.all <- c(z.all, z.others)
    }
    newz.grid <- expand.grid(z.all)
    z1save <- newz.grid[, 1]
    z2save <- newz.grid[, 2]
    colnames(newz.grid) <- colnames(Z)[ord]
    newz.grid <- newz.grid[,colnames(Z)]

    if(!is.null(min.plot.dist)) {
        mindists <- rep(NA, nrow(newz.grid))
        for(k in seq_along(mindists)) {
            pt <- as.numeric(newz.grid[k,c(colnames(Z)[ord[1]],colnames(Z)[ord[2]])])
            dists <- fields::rdist(matrix(pt, nrow = 1), Z[, c(colnames(Z)[ord[1]],colnames(Z)[ord[2]])])
            mindists[k] <- min(dists)
        }
    }

    if(preds.method == "approx") {
        preds <- ComputePostmeanHnew(fit = fit, y = y, Z = Z, X = X, Znew = newz.grid, sel = sel)
        preds.plot <- preds$postmean
        se.plot <- sqrt(diag(preds$postvar))
    } else if(preds.method == "samp") {
        stop("not yet implemented")
#         preds <- SampleHnew(Znew = newz.grid, fit = fit, Z = Z, X = X, y = y, ...)
#         preds.plot <- colMeans(preds)
#         se.plot <- apply(preds, 2, sd)
    }
    if(center) preds.plot <- preds.plot - mean(preds.plot)
    if(!is.null(min.plot.dist)) {
        preds.plot[mindists > min.plot.dist] <- NA
        se.plot[mindists > min.plot.dist] <- NA
    }
#     hgrid <- matrix(preds.plot, ngrid, ngrid, dimnames=list(z1=round(z1,2), z2=round(z2,2)))
#     se.grid <- matrix(se.plot, ngrid, ngrid, dimnames=list(z1=round(z1,2), z2=round(z2,2)))

    res <- dplyr::data_frame(z1 = z1save, z2 = z2save, est = preds.plot, se = se.plot)
}

#' Predict the exposure-response function at a new grid of points
#'
#' Predict the exposure-response function at a new grid of points
#'
#' @inheritParams kmbayes
#' @inheritParams ExtractEsts
#' @inheritParams SingVarRiskSummaries
#' @inheritParams PredictorResponseUnivar
#' @param z.pairs data frame showing which pairs of pollutants to plot
#' @param ngrid number of grid points in each dimension
#' @param verbose TRUE or FALSE: flag of whether to print intermediate output to the screen
#' @export
PredictorResponseBivar <- function(fit, y = NULL, Z = NULL, X = NULL, z.pairs = NULL, preds.method = "approx", ngrid = 50, q.fixed = 0.5, sel = NULL, min.plot.dist = 0.5, center = TRUE, z.names = colnames(Z), verbose = TRUE, ...) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
  }
  
  if (is.null(z.pairs)) {
    z.pairs <- expand.grid(z1 = 1:ncol(Z), z2 = 1:ncol(Z))
    z.pairs <- z.pairs[z.pairs$z1 < z.pairs$z2, ]
  }
  
  df <- dplyr::data_frame()
  for(i in 1:nrow(z.pairs)) {
    compute <- TRUE
    whichz1 <- z.pairs[i, 1] %>% unlist %>% unname
    whichz2 <- z.pairs[i, 2] %>% unlist %>% unname
    if(whichz1 == whichz2) compute <- FALSE
    z.name1 <- z.names[whichz1]
    z.name2 <- z.names[whichz2]
    names.pair <- c(z.name1, z.name2)
    if(nrow(df) > 0) { ## determine whether the current pair of variables has already been done
      completed.pairs <- df %>%
        dplyr::select_('variable1', 'variable2') %>%
        dplyr::distinct() %>%
        dplyr::transmute(z.pair = paste('variable1', 'variable2', sep = ":")) %>%
        unlist %>% unname
      if(paste(names.pair, collapse = ":") %in% completed.pairs | paste(rev(names.pair), collapse = ":") %in% completed.pairs) compute <- FALSE
    }
    if(compute) {
      if(verbose) message("Pair ", i, " out of ", nrow(z.pairs))
      res <- PredictorResponseBivarPair(fit = fit, y = y, Z = Z, X = X, whichz1 = whichz1, whichz2 = whichz2, preds.method = preds.method, ngrid = ngrid, q.fixed = q.fixed, sel = sel, min.plot.dist = min.plot.dist, center = center, z.names = z.names, ...)
      df0 <- res
      df0$variable1 <- z.name1
      df0$variable2 <- z.name2
      df0 %<>%
        dplyr::select_(~variable1, ~variable2, ~z1, ~z2, ~est, ~se)
      df <- dplyr::bind_rows(df, df0)
    }
  }
  df$variable1 <- factor(df$variable1, levels = z.names)
  df$variable2 <- factor(df$variable2, levels = z.names)
  df
}

#' Plot cross-sections of the bivariate predictor-response function
#' 
#' Function to plot the \code{h} function of a particular variable at different levels (quantiles) of a second variable
#' 
#' @export
#' @inheritParams kmbayes
#' @param pred.resp.df object obtained from running the function \code{\link{PredictorResponseBivar}}
#' @param qs vector of quantiles at which to fix the second variable
#' @param both_pairs flag indicating whether, if \code{h(z1)} is being plotted for z2 fixed at different levels, that they should be plotted in the reverse order as well (for \code{h(z2)} at different levels of z1) 
PredictorResponseBivarLevels <- function(pred.resp.df, Z = NULL, qs = c(0.25, 0.5, 0.75), both_pairs = TRUE) {
  var.pairs <- dplyr::distinct(dplyr::select_(pred.resp.df, ~variable1, ~variable2))
  if (both_pairs) {
    var.pairs.rev <- dplyr::data_frame(
      variable1 = var.pairs$variable2,
      
      variable2 = var.pairs$variable1
    )
    var.pairs <- rbind(var.pairs, var.pairs.rev)
  }
  
  df <- data.frame()
  for (i in 1:nrow(var.pairs)) {
    var1 <- as.character(unlist(var.pairs[i, "variable1"]))
    var2 <- as.character(unlist(var.pairs[i, "variable2"]))
    preds <- pred.resp.df[pred.resp.df$variable1 == var1 & pred.resp.df$variable2 == var2, ]
    if (nrow(preds) == 0) {
      preds <- pred.resp.df[pred.resp.df$variable1 == var2 & pred.resp.df$variable2 == var1, ]
    }
    
    ngrid <- sqrt(nrow(preds))
    preds.plot <- preds$est
    se.plot <- preds$se
    
    hgrid <- matrix(preds.plot, ngrid, ngrid)
    se.grid <- matrix(se.plot, ngrid, ngrid)
    z1 <- preds$z1[1:ngrid]
    z2 <- preds$z2[seq(1, by = ngrid, length.out = ngrid)]
    
    quants <- quantile(Z[, var2], qs)
    
    ## relation of z1 with outcome at different levels of z2
    se.grid.sub <- hgrid.sub <- matrix(NA, ngrid, length(qs))
    for (k in seq_along(quants)) {
      sub.sel <- which.min(abs(z2 - quants[k]))
      hgrid.sub[, k] <- hgrid[, sub.sel]
      se.grid.sub[, k] <- se.grid[, sub.sel]
    }
    colnames(hgrid.sub) <- colnames(se.grid.sub) <- paste0("q", seq_along(qs))
    hgrid.df <- tidyr::gather(data.frame(hgrid.sub), quantile, 'est', convert = TRUE)
    se.grid.df <- tidyr::gather(data.frame(se.grid.sub), quantile, 'se')
    
    df.curr <- data.frame(variable1 = var1, variable2 = var2, z1 = z1, quantile = factor(hgrid.df$quantile, labels = qs), est = hgrid.df$est, se = se.grid.df$se, stringsAsFactors = FALSE)
    df <- rbind(df, df.curr)
  }
  df <- dplyr::tbl_df(df) %>%
    dplyr::arrange_(~variable1, ~variable2)
  df
}
