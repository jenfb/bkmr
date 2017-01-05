devtools::load_all()

## gaussian

set.seed(111)
dat <- SimData(n = 200, M = 4)

set.seed(111)
orig <- kmbayes(y = dat$y, Z = dat$Z, X = dat$X, iter = 5000, verbose = FALSE, varsel = TRUE, est.h = TRUE)
set.seed(111)
new <- kmbayes(y = dat$y, Z = dat$Z, X = dat$X, iter = 5000, verbose = FALSE, varsel = TRUE, est.h = FALSE)

with(orig, difftime(time2, time1))
with(new, difftime(time2, time1))

## binomial

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
Xnew <- rbind(0)

tmp <- SamplePred(fit = new, Znew = Znew, Xnew = Xnew)
head(tmp)