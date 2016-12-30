devtools::load_all()
library(magrittr)

## generate data ####

seed <- 1234
#seed <- 8
#seed <- 31
set.seed(seed)
n = 200
M = 5
sigsq.true = 1
beta.true = 0.1
Z <- matrix(runif(n * M, -1, 1), n, M)
X <- cbind(3*cos(Z[, 1]) + 2*rnorm(n))
eps <- rnorm(n, sd = sqrt(sigsq.true))
#h <- 4*plogis(Z[, 1], 0, 0.3)
h <- (2*Z[, 1] + 0.5) ^ 2
plot(Z[, 1], h)

mu <- X * beta.true + h
y <- drop(mu + eps)
  ystar <- y
  y <- ifelse(ystar > 0, 1, 0)
dat <- list(n = n, M = M, sigsq.true = sigsq.true, beta.true = beta.true, Z = Z, h = h, X = X, y = y)
  dat$ystar <- ystar
datp <- dat
#datp <- SimData(hfun = 1, family = "binomial", beta.true = 0.1)
y <- datp$y
Z <- datp$Z
X <- datp$X
mu <- datp %$% drop(h + X*beta.true)
#str(datp)
table(y)
z1ord <- order(Z[, 1])

plot(Z[, 1], datp$ystar)
abline(lm(datp$h ~ Z[, 1]), col = "blue")
lines(Z[z1ord, 1], predict(lm(datp$h ~ splines::ns(Z[, 1], 3)))[z1ord], col = "blue")
lines(Z[z1ord, 1], datp$h[z1ord], col = "red", lwd = 2)

## linear probit model ####

probitfit0 <- glm(y ~ Z + X, family = binomial(link = "probit"))

tmp <- predict(probitfit0)
plot(tmp, datp$ystar)
abline(0, 1, col = "red")

hpred <- drop(coef(probitfit0)["(Intercept)"] + Z %*% coef(probitfit0)[grep("Z", names(coef(probitfit0)))])
plot(hpred, datp$h)
abline(0, 1, col = "red")

## BKM probit model ####

# fitpr_gam0 <- kmbayes(iter = 5000, y = y, Z = Z, X = X, family = "binomial", varsel = TRUE, control.params = list(verbose_show_ests = TRUE, r.prior = "gamma"))
# fitpr <- fitpr_gam0
# 
# fitpr_gam1 <- kmbayes(iter = 5000, y = y, Z = Z, X = X, family = "binomial", varsel = TRUE, control.params = list(verbose_show_ests = TRUE, r.prior = "gamma", r.jump2 = 1))
# fitpr <- fitpr_gam1

fitpr_gam2 <- kmbayes(iter = 5000, y = y, Z = Z, X = X, family = "binomial", varsel = TRUE, control.params = list(verbose_show_ests = TRUE, r.prior = "gamma", r.jump2 = 2))
fitpr <- fitpr_gam2

fitpr_iu <- kmbayes(iter = 5000, y = y, Z = Z, X = X, family = "binomial", varsel = TRUE, control.params = list(verbose_show_ests = TRUE, r.prior = "invunif"))
fitpr <- fitpr_iu

fitpr_iu2 <- kmbayes(iter = 5000, y = y, Z = Z, X = X, family = "binomial", varsel = TRUE, control.params = list(verbose_show_ests = TRUE, r.prior = "invunif", r.jump1 = 3))
fitpr <- fitpr_iu2

## with knots matrix

sel <- with(fitpr, seq(floor(iter/2) + 1, iter))
ests <- ExtractEsts(fitpr)
#summary(fitpr)

TracePlot(fit = fitpr, par = "beta")
TracePlot(fit = fitpr, par = "sigsq.eps")
TracePlot(fit = fitpr, par = "r", comp = 1)
TracePlot(fit = fitpr, par = "r", comp = 2)
TracePlot(fit = fitpr, par = "r", comp = 3)
TracePlot(fit = fitpr, par = "r", comp = 4)
TracePlot(fit = fitpr, par = "r", comp = 5)
TracePlot(fit = fitpr, par = "h", comp = 20)
TracePlot(fit = fitpr, par = "ystar", comp = 1)

## beta
cbind(truth = datp$beta.true, ests$beta) %>% t()

## h
hhat <- ests$h[, "mean"]
plot(hhat, datp$h)
abline(0, 1, col = "red")
#plot(Z[z1ord, 1], (datp$ystar - datp$X*0.1)[z1ord], col = "blue")
plot(Z[z1ord, 1], hhat[z1ord], ylim = range(c(hhat, datp$h)), col = ifelse(y[z1ord] == 1, "green", "blue"))
lines(Z[z1ord, 1], datp$h[z1ord], col = "red")

## ystar
ystar_hat <- colMeans(fitpr$ystar[sel, ])
plot(ystar_hat, datp$ystar)
abline(0, 1, col = "red")
par(mfrow = c(1,2))
plot(datp$y, datp$ystar)
plot(datp$y, ystar_hat)
par(mfrow = c(1,1))
par(mfrow = c(1,2))
hist(datp$ystar)
hist(ystar_hat)
par(mfrow = c(1,1))
par(mfrow = c(1,2))
plot(Z[z1ord, 1], datp$ystar[z1ord], col = ifelse(y[z1ord] == 1, "green", "blue"))
lines(Z[z1ord, 1], datp$h[z1ord], col = "red")
plot(Z[z1ord, 1], ests$ystar[z1ord], col = ifelse(y[z1ord] == 1, "green", "blue"))
lines(Z[z1ord, 1], hhat[z1ord], col = "red")
par(mfrow = c(1,1))

## phat
phat1 <- colMeans(pnorm(fitpr$h.hat[sel, ] + t(X %*% fitpr$beta[sel, ])))
phat2 <- colMeans(fitpr$ystar[sel, ] > 0)
table(phat2)

if (FALSE) {

  iter = 1000; family = "binomial"; id = NULL; verbose = TRUE; Znew = NULL; starting.values = list(); control.params = list(); varsel = FALSE; groups = NULL; knots = NULL; ztest = NULL; rmethod = "varying"
 
  s <- 2
  beta = chain$beta[s-1,]
  Vinv = Vcomps$Vinv
  Xbeta <-  drop(X %*% beta)
  lower <- ifelse(y == 1, 0, -Inf)
  upper <- ifelse(y == 0, 0,  Inf)
  time <- system.time(samp <- tmvtnorm::rtmvnorm(1, mean = Xbeta, H = Vinv, lower = lower, upper = upper, algorithm = "gibbs", start.value = chain$ystar[s - 1, ]))
  time
  head(drop(samp))
  
}