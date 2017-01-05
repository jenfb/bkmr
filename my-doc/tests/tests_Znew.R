devtools::load_all()

set.seed(111)
dat <- SimData(n = 50, M = 4)
y <- dat$y
Z <- dat$Z
X <- dat$X

Znew <- rbind(c(0,0,0,0), c(1,1,1,1))

set.seed(111)
fit0 <- kmbayes(y = y, Z = Z, X = X, iter = 1000, verbose = FALSE, varsel = TRUE, Znew = Znew)

ExtractEsts(fit0)$hnew

ComputePostmeanHnew(fit0, Znew = Znew)

fit0

summary(fit0)
