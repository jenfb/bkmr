CalcR2 <- function(fit, covar) {
    ests <- ExtractEsts(fit)

    ## http://en.wikipedia.org/wiki/Coefficient_of_determination
    preds <- with(ests, h[, "mean"] + covar %*% beta[, "mean"])
    SStot <- sum((y - mean(y))^2)
    SSres <- sum((y - preds)^2)
    R2 <- 1 - SSres/SStot
    R2
}
