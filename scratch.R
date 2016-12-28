devtools::load_all()

set.seed(1234)
datp <- SimData(hfun = 1, family = "binomial", beta.true = 1)
y <- datp$y
Z <- datp$Z
X <- datp$X
#str(datp)
table(y)

plot(Z[, 1], datp$ystar)

probitfit0 <- glm(y ~ Z + X, family = binomial(link = "probit"))
tmp <- predict(probitfit0)
plot(tmp, datp$ystar)
abline(0, 1, col = "red")

fitpr <- kmbayes(iter = 10000, y = y, Z = Z, X = X, family = "binomial", varsel = TRUE, control.params = list(verbose_show_ests = TRUE))

summary(fitpr)

hhat <- ExtractEsts(fitpr)$h[, "mean"]
plot(hhat, datp$h)
abline(0, 1, col = "red")

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