devtools::load_all()

## generate data ####

#seed <- 1234
#seed <- 8
seed <- 31
set.seed(seed)
datp <- SimData(hfun = 1, family = "binomial", beta.true = 0.1)
y <- datp$y
Z <- datp$Z
X <- datp$X
mu <- datp %$% drop(h + X*beta.true)
y
#str(datp)
table(y)
z1ord <- order(Z[, 1])

plot(Z[, 1], datp$ystar)
abline(lm(datp$h ~ Z[, 1]), col = "blue")
lines(Z[z1ord, 1], predict(lm(datp$h ~ splines::ns(Z[, 1], 3)))[ord], col = "blue")
lines(Z[z1ord, 1], datp$h[ord], col = "red", lwd = 2)

## linear probit model ####

probitfit0 <- glm(y ~ Z + X, family = binomial(link = "probit"))

tmp <- predict(probitfit0)
plot(tmp, datp$ystar)
abline(0, 1, col = "red")

hpred <- drop(coef(probitfit0)["(Intercept)"] + Z %*% coef(probitfit0)[grep("Z", names(coef(probitfit0)))])
plot(hpred, datp$h)
abline(0, 1, col = "red")

## BKM probit model ####

#fitpr <- kmbayes(iter = 10000, y = y, Z = Z, X = X, family = "binomial", varsel = TRUE, control.params = list(verbose_show_ests = TRUE, r.prior = "gamma"))
#fitpr <- kmbayes(iter = 5000, y = y, Z = Z, X = X, family = "binomial", varsel = TRUE, control.params = list(verbose_show_ests = TRUE, r.prior = "invunif"))
sel <- with(fitpr, seq(floor(iter/2) + 1, iter))
ests <- ExtractEsts(fitpr)
#summary(fitpr)

## beta
cbind(truth = datp$beta.true, ests$beta) %>% t()

## h
hhat <- ests$h[, "mean"]
plot(hhat, datp$h)
abline(0, 1, col = "red")
#plot(Z[z1ord, 1], (datp$ystar - datp$X*0.1)[z1ord], col = "blue")
plot(Z[z1ord, 1], hhat[z1ord], ylim = range(c(hhat, datp$h)))
lines(Z[z1ord, 1], datp$h[z1ord], col = "red")

## ystar
ystar_hat <- colMeans(fitpr$ystar[sel, ])
plot(ystar_hat, datp$ystar)
abline(0, 1, col = "red")

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