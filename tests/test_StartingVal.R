##context("starting values validation")

library(bkmr)
library(testthat)
set.seed(111)
dat <- SimData(n = 50, M = 4)
y <- dat$y
Z <- dat$Z
X <- dat$X
sval <- list(h.hat = 1, beta = 1, sigsq.eps = 1, r = 1, lambda = 10, delta = 1)
sval1 <- list(h.hat = seq(1:50), beta = 1, sigsq.eps = 1, r = 1, lambda = 10, delta = 1)
X2<-cbind(X, seq(1:50))
test_that("starting.values validation-length of list elements", {
  expect_message (testfit <- kmbayes(y = y, Z = Z, X = X, iter = 100, starting.values = sval), 
          'h.hat should be a vector of length equal to number of rows in Y.  Input will be repeated or truncated as necessary.', fixed=TRUE)
  expect_message (testfit <- kmbayes(y = y, Z = Z, X = X, iter = 100, starting.values = sval1), 
            'delta should be a vector of length equal to the number of columns of Z.  Input will be repeated or truncated as necessary.', fixed=TRUE)
  expect_message (testfit <- kmbayes(y = y, Z = Z, X = X2, iter = 100, starting.values = sval), 
              'beta should be a vector of length equal to the number of columns of X.  Input will be repeated or truncated as necessary.', fixed=TRUE)
})
sval$h.hat<-seq(1:50)
sval$delta<-rep(10,4)
sval1$h.hat<-rep(0, 50)
sval2 <- list(h.hat = seq(1:50), beta = 1, sigsq.eps = 1, r = 1, lambda = 10, delta = 1)
sval2$sigsq.eps <- 0
sval3 <- list(h.hat = seq(1:50), beta = 1, sigsq.eps = 1, r = 1, lambda = 10, delta = 1)
sval3$lambda <- 0
test_that("check individual values of vectors", {
  expect_error (testfit <- kmbayes(y = y, Z = Z, X = X, iter = 100, starting.values = sval), 
                'starting.values$delta == 1 || starting.values$delta == 0 is not TRUE', fixed=TRUE)
  expect_error (testfit <- kmbayes(y = y, Z = Z, X = X, iter = 100, starting.values = sval1), 
                'starting.values$h.hat > 0 are not all TRUE', fixed=TRUE)
  expect_error (testfit <- kmbayes(y = y, Z = Z, X = X, iter = 100, starting.values = sval2), 
                'starting.values$sigsq.eps > 0 is not TRUE', fixed=TRUE)
  expect_error (testfit <- kmbayes(y = y, Z = Z, X = X, iter = 100, starting.values = sval3), 
                'starting.values$lambda > 0 is not TRUE', fixed=TRUE)
})
sval5 <- list(h.hat = seq(1:50), beta = 1, sigsq.eps = 1, r = 1, lambda = 10, delta = rep(0,4))
sval6 <- list(h.hat = seq(1:50), beta = 1, sigsq.eps = 1, r = 0, lambda = 10, delta = rep(0,4))
sval7 <- list(h.hat = seq(1:50), beta = 1, sigsq.eps = 1, r = seq(1:4), lambda = 10, delta = rep(0,4))
sval8 <- list(h.hat = seq(1:50), beta = 1, sigsq.eps = 1, r = seq(1:4), lambda = 10, delta = rep(0,4))
sval8$r <- rep(0,4)
test_that("r, with and without varsel", {
  expect_message (testfit <- kmbayes(y = y, Z = Z, X = X, iter = 100, starting.values = sval5, varsel=TRUE), 
                'r should be a vector of length equal to the number of columns of Z.  Input will be repeated or truncated as necessary.', fixed=TRUE)
  expect_error (testfit <- kmbayes(y = y, Z = Z, X = X, iter = 100, starting.values = sval6, varsel=FALSE), 
                'starting.values$r > 0 is not TRUE', fixed=TRUE)
  expect_message (testfit <- kmbayes(y = y, Z = Z, X = X, iter = 100, starting.values = sval7, varsel=FALSE), 
                'r should a scalar.  Vector input will be truncated.', fixed=TRUE)
  expect_error (testfit <- kmbayes(y = y, Z = Z, X = X, iter = 100, starting.values = sval8, varsel=TRUE), 
                'starting.values$r > 0 are not all TRUE', fixed=TRUE)

})



