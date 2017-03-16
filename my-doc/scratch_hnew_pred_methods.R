set.seed(111)
dat <- SimData(n = 500, M = 4)
y <- dat$y
Z <- dat$Z
X <- dat$X

Xnew <- matrix(0, nrow(Z), ncol(X))

set.seed(111)
fit0 <- kmbayes(y = y, Z = Z, X = X, iter = 2000, verbose = FALSE, varsel = TRUE)

ests_full <- ExtractEsts(fit0)$h[, c("mean", "sd")]

## add Vcomps
sel <- with(fit0, seq(floor(iter/2) + 1, iter, 10))
sel <- unique(floor(sel))
s0 <- system.time(
  Vinv <- with(fit0, lapply(sel, function(s) makeVcomps(r = r[s, ], lambda = lambda[s, ], Z, data.comps)$Vinv))
)
attr(Vinv, "sel") <- sel
fit0V <- fit0
fit0V$Vinv <- Vinv
print(object.size(fit0), units = "Mb")
print(object.size(fit0V), units = "Mb")

s1 <- system.time(
  samps <- SamplePred(fit0, Xnew = Xnew)
)

ests_samp <- samps %>%
  {cbind(mean = colMeans(.), sd = apply(., 2, sd))}

s2 <- system.time(
  ests_approx <- ComputePostmeanHnew(fit0) %$%
    cbind(mean = postmean, sd = sqrt(diag(postvar)))
)

s3 <- system.time(
  ests_approx2 <- ComputePostmeanHnew2(fit0) %$%
    cbind(mean = postmean, sd = sqrt(diag(postvar)))
)

s4 <- system.time(
  ests_approx2 <- ComputePostmeanHnew2(fit0V) %$%
    cbind(mean = postmean, sd = sqrt(diag(postvar)))
)

s1
s2
s3
s4

head(ests_full)
head(ests_samp)
head(ests_approx)
head(ests_approx2)

par(mfrow = c(1, 2))
plot(ests_samp[, "mean"], ests_approx2[, "mean"])
abline(0, 1, col = "red")
plot(ests_samp[, "sd"], ests_approx2[, "sd"])
abline(0, 1, col = "red")

par(mfrow = c(1, 2))
plot(ests_samp[, "mean"], ests_approx[, "mean"])
abline(0, 1, col = "red")
plot(ests_samp[, "sd"], ests_approx[, "sd"])
abline(0, 1, col = "red")

par(mfrow = c(1, 2))
plot(ests_approx[, "mean"], ests_approx2[, "mean"])
abline(0, 1, col = "red")
plot(ests_approx[, "sd"], ests_approx2[, "sd"])
abline(0, 1, col = "red")

