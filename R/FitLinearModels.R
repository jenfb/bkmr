#' @export
SinglePolEsts <- function(y, expos, covar) {
    coefs <- matrix(NA, ncol(expos), 4)
    for(i in 1:ncol(expos)) {
        x <- expos[, i] %>% unlist %>% unname
        lm1 <- lm(y ~ x + covar)
        coefs[i, ] <- summary(lm1)$coef["x",]
    }
    df <- data_frame(exposure = colnames(expos),
                     est = coefs[, 1],
                     se = coefs[, 2],
                     tstat = coefs[, 3],
                     p = coefs[, 4])
    df
}

#' @export
MultiPolEsts <- function(lmall, y, expos, covar) {
    if (is.null(lmall)) {
        lmall <- lm(y ~ expos + covar)
    }
    sel <- grep("expos", names(coef(lmall)), value = TRUE)
    coefs <- summary(lmall)$coef[sel, ]
    df <- data_frame(exposure = colnames(expos),
                     est = coefs[, 1],
                     se = coefs[, 2],
                     tstat = coefs[, 3],
                     p = coefs[, 4])
    df
}
