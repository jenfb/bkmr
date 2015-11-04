context("argument validation")

test_yvalidation <- function() {
  library(bkmr)
  expect_error(kmbayes(y = vector('numeric'), Z = Z, X = X, iter = 100), 
                'length(y) > 0 is not TRUE', fixed=TRUE)
  expect_error(kmbayes(y = c(1, 'hello'), Z = Z, X = X, iter = 100), 
               'is.numeric(y) is not TRUE', fixed=TRUE)
  expect_error(kmbayes(y = c(1,NA,0), Z = Z, X = X, iter = 100), 
               'anyNA(y) == FALSE is not TRUE', fixed=TRUE)
}

test_Zvalidation <- function() {
  library(bkmr)
  ##expect_error(kmbayes(y = y, Z = vector('numeric'), X = X, iter = 100), 'ncol(Z) > 0 is not TRUE', fixed=TRUE)
  expect_error(kmbayes(y = y, Z = c(1, 'hello'), X = X, iter = 100), 
               'is.numeric(Z) is not TRUE', fixed=TRUE)
  expect_error(kmbayes(y = y, Z = cbind(matrix(c(1,1), 2), c(1,1)), X = X, iter = 100), 
               'nrow(Z) == length(y) is not TRUE', fixed=TRUE)
  expect_error(kmbayes(y = c(1,1), Z = cbind(matrix(c(1,NA), 2), c(1,1)), X = X, iter = 100), 
              'anyNA(Z) == FALSE is not TRUE', fixed=TRUE)
}

test_Xvalidation <- function() {
  library(bkmr)
  ##expect_error(kmbayes(y = y, Z = Z, X = vector('numeric'), iter = 100),  'ncol(X) > 0 is not TRUE', fixed=TRUE)
  expect_error(kmbayes(y = y, Z = Z, X = c(1, 'hello'), iter = 100), 
               'is.numeric(X) is not TRUE', fixed=TRUE)
  expect_error(kmbayes(y = y, Z = Z, X = cbind(matrix(c(1,1), 2), c(1,1)), iter = 100), 
               'nrow(X) == length(y) is not TRUE', fixed=TRUE)
  expect_error(kmbayes(y = c(1,1), Z = c(1,1), X = cbind(matrix(c(1,NA), 2), c(1,1)), iter = 100), 
               'anyNA(X) == FALSE is not TRUE', fixed=TRUE)
  ## create a fail
 ## expect_error(kmbayes(y = c(1,1), Z = c(1,1), X = cbind(matrix(c(1,NA), 2), c(1,1)), iter = 100), 
   ##            'hello!', fixed=TRUE)
 
}

test_defaultreset <- function() {
  library(bkmr)
  set.seed(111)
  dat <- SimData(n = 50, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X
  expect_message(kmbayes(y=y, X=X, Z=Z, iter=0), 
                 'invalid input for iter, resetting to default value 1000', fixed=TRUE)
  expect_message(kmbayes(y=y, X=X, Z=Z, iter=10, family= 'wrong'), 
                 'not yet implemented, resetting family to default gaussian', fixed=TRUE)
  expect_message(kmbayes(y=y, X=X, Z=Z, iter=10, verbose = 'wrong'), 
                 'invalid value for verbose, resetting to default FALSE', fixed=TRUE)
  expect_message(kmbayes(y=y, X=X, Z=Z, iter=10, varsel = 'wrong'), 
                 'invalid value for varsel, resetting to default FALSE', fixed=TRUE)
}

test_argsRemaining <- function() {
  library(bkmr)
  set.seed(111)
  dat <- SimData(n = 50, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X
 ##id length y, id not missing, id not knots
  expect_error(kmbayes(y=head(y, 10), X=head(X, 10), Z=head(Z, 10), id=c(1,1)), 
                 'length(id) == length(y) is not TRUE', fixed=TRUE)
  expect_error(kmbayes(y=c(1,1), X=X, Z=c(1,1), id=c(1,NA)), 
                 '', fixed=TRUE)
  
}


