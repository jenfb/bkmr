devtools::load_all()
library(magrittr)

## Gaussian outcome

set.seed(111)
dat <- SimData(n = 50, M = 4)
y <- dat$y
Z <- dat$Z
X <- dat$X

Znew <- rbind(c(0,0,0,0), c(1,1,1,1))

set.seed(111)
fit0 <- kmbayes(y = y, Z = Z, X = X, iter = 10000, verbose = FALSE, varsel = TRUE, Znew = Znew)

ests_samp <- ExtractEsts(fit0)$hnew
ests_samp2 <- SamplePred(fit0, Znew = Znew, Xnew = matrix(0, nrow(Znew), ncol(datp$X))) %>%
{cbind(colMeans(.), apply(., 2, sd))}
system.time(
  ests_approx <- ComputePostmeanHnew(fit0, Znew = Znew) %$%
    cbind(mean = postmean, sd = sqrt(diag(postvar)))
)
system.time(
  ests_approx2 <- ComputePostmeanHnew2(fit0, Znew = Znew) %$%
    cbind(mean = postmean, sd = sqrt(diag(postvar)))
)

ests_samp
ests_samp2
ests_approx
ests_approx2

fit0

summary(fit0)

## Binomial outcome

set.seed(123)
n <- 200 ## number of observations
M <- 4 ## number of exposure variables
beta.true <- 0.1
Z <- matrix(runif(n * M, -1, 1), n, M)
x <- 3*cos(Z[, 1]) + 2*rnorm(n)
hfun <- function(z) (2*z + 0.5) ^ 2
h <- hfun(Z[, 1]) ## only depends on z1
## generate using latent normal representation
eps <- rnorm(n)
ystar <- x * beta.true + h + eps
y <- ifelse(ystar > 0, 1, 0)
datp <- list(n = n, M = M, beta.true = beta.true, Z = Z, h = h, X = cbind(x), y = y, eps = eps, ystar = ystar)
rm(n, M, beta.true, Z, x, h, eps, y, ystar)

Znew <- rbind(c(0,0,0,0), c(1,1,1,1))

set.seed(123)
fit0 <- kmbayes(y = datp$y, Z = datp$Z, X = datp$X, iter = 5000, verbose = TRUE, varsel = TRUE, Znew = Znew, family = "binomial", control.params = list(r.jump2 = 0.5))

ests_samp <- ExtractEsts(fit0)$hnew
system.time(
  ests_samp2 <- SamplePred(fit0, Znew = Znew, Xnew = matrix(0, nrow(Znew), ncol(datp$X))) %>%
    {cbind(colMeans(.), apply(., 2, sd))}
)
system.time(
  ests_approx <- ComputePostmeanHnew(fit0, Znew = Znew) %$%
    cbind(mean = postmean, sd = sqrt(diag(postvar)))
)
system.time(
  ests_approx2 <- ComputePostmeanHnew2(fit0, Znew = Znew) %$%
    cbind(mean = postmean, sd = sqrt(diag(postvar)))
)

ests_samp
ests_samp2
ests_approx
ests_approx2

fit0

summary(fit0)
