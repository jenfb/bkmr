#' @export
ExposureResponseUnivarPol <- function(whichpol = 1, fit, y, expos, covar, preds.method = "approx", ngrid = 50, q.fixed = 0.5, min.plot.dist = Inf, center = TRUE, expos.names = colnames(expos), ...) {

    if(ncol(expos) < 2) stop("requires there to be at least 2 exposure variables")

    if(is.null(expos.names)) {
        colnames(expos) <- paste0("espos", 1:ncol(expos))
    } else {
        colnames(expos) <- expos.names
    }

    ord <- c(whichpol, setdiff(1:ncol(expos), whichpol))
    z1 <- seq(min(expos[,ord[1]]), max(expos[,ord[1]]), length = ngrid)
    z.others <- lapply(2:ncol(expos), function(x) quantile(expos[,ord[x]], q.fixed))
    z.all <- c(list(z1), z.others)
    newz.grid <- expand.grid(z.all)
    colnames(newz.grid) <- colnames(expos)[ord]
    newz.grid <- newz.grid[,colnames(expos)]

    if(!is.null(min.plot.dist)) {
        mindists <- rep(NA,nrow(newz.grid))
        for(i in seq_along(mindists)) {
            pt <- as.numeric(newz.grid[i, colnames(expos)[ord[1]]])
            dists <- fields::rdist(matrix(pt, nrow = 1), expos[, colnames(expos)[ord[1]]])
            mindists[i] <- min(dists)
        }
    }

    if(preds.method == "approx") {
        preds <- ComputePostmeanHnew(fit = fit, y = y, expos = expos, covar = covar, exposNew = newz.grid)
        preds.plot <- preds$postmean
        se.plot <- sqrt(diag(preds$postvar))
    } else if(preds.method == "samp") {
        stop("not yet implemented")
#         preds <- SampleHnew(Znew = newz.grid, fit = fit, expos = expos, covar = covar, y = y, ...)
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

#' @export
ExposureResponseUnivar <- function(fit, y, expos, covar, pollutants = 1:ncol(expos), preds.method = "approx", ngrid = 50, q.fixed = 0.5, min.plot.dist = Inf, center = TRUE, expos.names = colnames(expos), ...) {
    df <- dplyr::data_frame()
    for(i in pollutants) {
        res <- ExposureResponseUnivarPol(whichpol = i, fit = fit, y = y, expos = expos, covar = covar, preds.method = preds.method, ngrid = ngrid, q.fixed = q.fixed, min.plot.dist = min.plot.dist, center = center, expos.names = expos.names, ...)
        df0 <- dplyr::mutate(res, exposure = expos.names[i]) %>%
            dplyr::select(exposure, z, est, se)
        df <- dplyr::bind_rows(df, df0)
    }
    df <- dplyr::mutate(df, exposure = factor(exposure, levels = expos.names[pollutants]))
}




#' @export
ExposureResponseBivarPair <- function(fit, y, expos, covar, whichpol1 = 1, whichpol2 = 2, whichpol3, preds.method = "approx", prob = 0.5, q.fixed = 0.5, ngrid = 50, min.plot.dist = 0.5, center = TRUE, ...) {
    if(ncol(expos) < 3) stop("requires there to be at least 3 exposure variables")

    if(is.null(colnames(expos))) colnames(expos) <- paste0("espos", 1:ncol(expos))

    if(missing(whichpol3)) {
        ord <- c(whichpol1, whichpol2, setdiff(1:ncol(expos), c(whichpol1, whichpol2)))
    } else {
        ord <- c(whichpol1, whichpol2, whichpol3, setdiff(1:ncol(expos), c(whichpol1, whichpol2, whichpol3)))
    }
    z1 <- seq(min(expos[,ord[1]]), max(expos[,ord[1]]), length=ngrid)
    z2 <- seq(min(expos[,ord[2]]), max(expos[,ord[2]]), length=ngrid)
    z3 <- quantile(expos[, ord[3]], probs = prob)
    z.all <- c(list(z1), list(z2), list(z3))
    if(ncol(expos) > 3) {
        z.others <- lapply(4:ncol(expos), function(x) quantile(expos[,ord[x]], q.fixed))
        z.all <- c(z.all, z.others)
    }
    newz.grid <- expand.grid(z.all)
    z1save <- newz.grid[, 1]
    z2save <- newz.grid[, 2]
    colnames(newz.grid) <- colnames(expos)[ord]
    newz.grid <- newz.grid[,colnames(expos)]

    if(!is.null(min.plot.dist)) {
        mindists <- rep(NA, nrow(newz.grid))
        for(k in seq_along(mindists)) {
            pt <- as.numeric(newz.grid[k,c(colnames(expos)[ord[1]],colnames(expos)[ord[2]])])
            dists <- fields::rdist(matrix(pt, nrow = 1), expos[, c(colnames(expos)[ord[1]],colnames(expos)[ord[2]])])
            mindists[k] <- min(dists)
        }
    }

    if(preds.method == "approx") {
        preds <- ComputePostmeanHnew(fit = fit, y = y, expos = expos, covar = covar, exposNew = newz.grid)
        preds.plot <- preds$postmean
        se.plot <- sqrt(diag(preds$postvar))
    } else if(preds.method == "samp") {
        stop("not yet implemented")
#         preds <- SampleHnew(Znew = newz.grid, fit = fit, expos = expos, covar = covar, y = y, ...)
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

#' @export
ExposureResponseBivar <- function(fit, y, expos, covar, pollutant.pairs = subset(expand.grid(x = 1:ncol(expos), y = 1:ncol(expos)), x < y), preds.method = "approx", ngrid = 50, q.fixed = 0.5, min.plot.dist = 0.5, center = TRUE, expos.names = colnames(expos), quiet = FALSE, ...) {

    df <- dplyr::data_frame()
    for(i in 1:nrow(pollutant.pairs)) {
        compute <- TRUE
        whichpol1 <- pollutant.pairs[i, 1] %>% unlist %>% unname
        whichpol2 <- pollutant.pairs[i, 2] %>% unlist %>% unname
        if(whichpol1 == whichpol2) compute <- FALSE
        expos.name1 <- expos.names[whichpol1]
        expos.name2 <- expos.names[whichpol2]
        names.pair <- c(expos.name1, expos.name2)
        if(nrow(df) > 0) { ## determine whether the current pair of pollutants has already been done
            completed.pairs <- df %>%
                dplyr::select(exposure1, exposure2) %>%
                dplyr::distinct() %>%
                dplyr::transmute(expos.pair = paste(exposure1, exposure2, sep = ":")) %>%
                unlist %>% unname
            if(paste(names.pair, collapse = ":") %in% completed.pairs | paste(rev(names.pair), collapse = ":") %in% completed.pairs) compute <- FALSE
        }
        if(compute) {
            if(!quiet) message("Pair ", i, " out of ", nrow(pollutant.pairs))
            res <- ExposureResponseBivarPair(fit = fit, y = y, expos = expos, covar = covar, whichpol1 = whichpol1, whichpol2 = whichpol2, preds.method = preds.method, ngrid = ngrid, q.fixed = q.fixed, min.plot.dist = min.plot.dist, center = center, expos.names = expos.names, ...)
            df0 <- dplyr::mutate(res, exposure1 = expos.name1, exposure2 = expos.name2) %>%
                dplyr::select(exposure1, exposure2, z1, z2, est, se)
            df <- dplyr::bind_rows(df, df0)
        }
    }
    df <- dplyr::mutate(df, exposure1 = factor(exposure1, levels = expos.names), exposure2 = factor(exposure2, levels = expos.names))
}

#' Function to plot the exposure response function of a particular pollutant at different levels (quantiles) of a second pollutant
#' @export
#' @param expos.resp.df object obtained from running the function \code{ExposureResponseBivar()}
#' @param qs vector of quantiles of the second pollutant
ExposureResponseBivarLevels <- function(expos.resp.df, expos, qs = c(0.25, 0.5, 0.75)) {
    pol.pairs <- dplyr::distinct(dplyr::select(expos.resp.df, exposure1, exposure2))

    df <- data.frame()
    for (i in 1:nrow(pol.pairs)) {
        exp1 <- as.character(unlist(pol.pairs[i, "exposure1"]))
        exp2 <- as.character(unlist(pol.pairs[i, "exposure2"]))
        preds <- subset(expos.resp.df, exposure1 == exp1 & exposure2 == exp2)

        ngrid <- sqrt(nrow(preds))
        preds.plot <- preds$est
        se.plot <- preds$se

        hgrid <- matrix(preds.plot, ngrid, ngrid)
        se.grid <- matrix(se.plot, ngrid, ngrid)
        z1 <- preds$z1[1:ngrid]
        z2 <- preds$z2[seq(1, by = ngrid, length.out = ngrid)]

        for (j in 1:2) {
            if (j == 1) {
                quants <- quantile(expos[, exp2], qs)
            } else if (j == 2) { ## switch the roles of z1, z2
                exp1new <- exp2
                exp2 <- exp1
                exp1 <- exp1new

                z1new <- z2
                z2 <- z1
                z1 <- z1new

                quants <- quantile(expos[, exp1], qs)
                hgrid <- t(hgrid)
                se.grid <- t(se.grid)
            }

            ## relation of z1 with outcome at different levels of z2
            se.grid.sub <- hgrid.sub <- matrix(NA, ngrid, length(qs))
            for (k in seq_along(quants)) {
                sub.sel <- which.min(abs(z2 - quants[k]))
                hgrid.sub[, k] <- hgrid[, sub.sel]
                se.grid.sub[, k] <- se.grid[, sub.sel]
            }
            colnames(hgrid.sub) <- colnames(se.grid.sub) <- paste0("q", seq_along(qs))
            hgrid.df <- tidyr::gather(data.frame(hgrid.sub), quantile, est, convert = TRUE)
            se.grid.df <- tidyr::gather(data.frame(se.grid.sub), quantile, se)

            df.curr <- data.frame(exposure1 = exp1, exposure2 = exp2, z1 = z1, quantile = factor(hgrid.df$quantile, labels = qs), est = hgrid.df$est, se = se.grid.df$se, stringsAsFactors = FALSE)
            df <- rbind(df, df.curr)
        }
    }
    df <- dplyr::tbl_df(df) %>%
        dplyr::arrange(exposure1, exposure2)
    df
}
