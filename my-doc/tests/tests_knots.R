devtools::load_all()
library(ggplot2)

family <- "gaussian"
#family <- "binomial"

## example using a matrix of knots ####

n <- 500
M <- 5
sigsq.true <- ifelse(family == "gaussian", 0.05, 1)
beta.true <- 0.5
Z <- matrix(rnorm(n * M), n, M)
X <- cbind(3*cos(Z[, 1]) + 2*rnorm(n))
eps <- rnorm(n, sd = sqrt(sigsq.true))
h <- apply(Z, 1, function(z, ind = 1) 4*plogis(z[ind[1]], 0, 0.3))
eps <- rnorm(n, sd = sqrt(sigsq.true))
y <- drop(X*beta.true + h + eps)
if (family == "binomial") {
  ystar <- y
  y <- ifelse(ystar > 0, 1, 0)
}

des <- fields::cover.design(Z, nd = 100)
knots <- des$design

## compare timing to not using knots 
set.seed(111)
fit_full <- kmbayes(y = y, Z = Z, X = X, iter = 300, varsel = TRUE, family = family, verbose = FALSE)
with(fit_full, difftime(time2, time1))
fit_knots <- kmbayes(y = y, Z = Z, X = X, iter = 300, varsel = TRUE, family = family, verbose = FALSE, knots = knots)
with(fit_knots, difftime(time2, time1))

## with with a larger number of iterations
fit0 <- kmbayes(y = y, Z = Z, X = X, iter = 5000, varsel = TRUE, family = family, knots = knots, control.params = list(verbose_show_ests = TRUE))

summary(fit0)

TracePlot(fit = fit0, par = "beta")
ExtractPIPs(fit0)