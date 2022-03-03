CalcR2 <- function(fit, sel = NULL) {
    ests <- ExtractEsts(fit, sel = sel)
    y <- fit$y

    ## http://en.wikipedia.org/wiki/Coefficient_of_determination
    preds <- with(ests, h[, "mean"] + fit$X %*% beta[, "mean"])
    SStot <- sum((y - mean(y))^2)
    SSres <- sum((y - preds)^2)
    R2 <- 1 - SSres/SStot
    R2
}
