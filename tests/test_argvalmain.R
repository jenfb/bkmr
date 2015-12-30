##context("argument validation")

library(bkmr)
library(testthat)
test_that("y validation", {
  expect_error(testfit<-kmbayes(y = vector('numeric'), Z = Z, X = X, iter = 100), 
                'length(y) > 0 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = c(1, 'hello'), Z = Z, X = X, iter = 100), 
               'is.numeric(y) is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = c(1,NA,0), Z = Z, X = X, iter = 100), 
               'anyNA(y) == FALSE is not TRUE', fixed=TRUE)
})
set.seed(111)
dat <- SimData(n = 50, M = 4)
y <- dat$y
Z <- dat$Z
X <- dat$X
test_that ("X validation", {
  ##expect_error(testfit<-kmbayes(y = y, Z = Z, X = vector('numeric'), iter = 100),  'ncol(X) > 0 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = c(1, 'hello'), iter = 100), 
               'is.numeric(X) is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = cbind(matrix(c(1,1), 2), c(1,1)), iter = 100), 
               'nrow(X) == length(y) is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = c(1,1), Z = c(1,1), X = cbind(matrix(c(1,NA), 2), c(1,1)), iter = 100), 
               'anyNA(X) == FALSE is not TRUE', fixed=TRUE)
  ## create a fail
 ## expect_error(kmbayes(y = c(1,1), Z = c(1,1), X = cbind(matrix(c(1,NA), 2), c(1,1)), iter = 100), 
         ##      'hello!', fixed=TRUE)
  
})

test_that("reset invalid to default", {
  expect_message(testfit<-kmbayes(y=y, X=X, Z=Z, iter=0), 
                 'invalid input for iter, resetting to default value 1000', fixed=TRUE)
  expect_message(testfit<-kmbayes(y=y, X=X, Z=Z, iter=10, family= 'wrong'), 
                 'not yet implemented, resetting family to default gaussian', fixed=TRUE)
  expect_message(testfit<-kmbayes(y=y, X=X, Z=Z, iter=10, verbose = 'wrong'), 
                 'invalid value for verbose, resetting to default FALSE', fixed=TRUE)
  expect_message(testfit<-kmbayes(y=y, X=X, Z=Z, iter=10, varsel = 'wrong'), 
                 'invalid value for varsel, resetting to default FALSE', fixed=TRUE)
})

Zmiss<-Z
Zmiss[,1][Zmiss[,1] > 0] <- NA
myid<-c(seq(1, 48), c(1,2))
test_that("validate remaining arguments",  {
  ##id: length y, id not missing, id not knots
  expect_error(testfit<-kmbayes(y=head(y, 10), X=head(X, 10), Z=head(Z, 10), id=c(1,1)), 
               'length(id) == length(y) is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y=head(y, 10), X=head(X, 10), Z=head(Z, 10), id=c(1,1,1,2,NA,2,1,1,1,2)), 
               'anyNA(id) == FALSE is not TRUE', fixed=TRUE)
  ##this raises an actual error:
  ##expect_message(testfit<-kmbayes(y=y, X=X, Z=Z, id=myid, knots=head(Z, 10)), 
         ##        'knots cannot be specified with id, resetting knots to null', fixed=TRUE)
  ##Znew: numeric, same columns as Z, no missing
  expect_error(testfit<-kmbayes(y=head(y, 10), X=head(X, 10), Z=head(Z, 10), Znew=c("hello", 1)), 
               'is.numeric(Znew) is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y=head(y, 10), X=head(X, 10), Z=head(Z, 10), Znew=tail(y,10)), 
               'ncol(Znew) == ncol(Z) is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y=y, X=X, Z=Z, Znew=Zmiss), 
               'anyNA(Znew) == FALSE is not TRUE', fixed=TRUE)
  ##knots: numeric, same columns as Z, no missing
  expect_error(testfit<-kmbayes(y=head(y, 10), X=head(X, 10), Z=head(Z, 10), knots=c("hello", 1)), 
               'is.numeric(knots) is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y=head(y, 10), X=head(X, 10), Z=head(Z, 10), knots=tail(y,10)), 
               'ncol(knots) == ncol(Z) is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y=y, X=X, Z=Z, knots=Zmiss), 
               'anyNA(knots) == FALSE is not TRUE', fixed=TRUE)
  ##groups: only with varsel, numeric, same columns as Z, no missing
  expect_message(testfit<-kmbayes(y=y, X=X, Z=Z, groups=c(1,2,2,1)), 
                 'groups should only be specified if varsel=TRUE, resetting varsel to TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y=head(y, 10), X=head(X, 10), Z=head(Z, 10), varsel=TRUE, groups=c(1,1,1)),
               'length(groups) == ncol(Z) is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y=head(y, 10), X=head(X, 10), Z=head(Z, 10), varsel=TRUE, groups=c("hello",1,1,1)), 
               'is.numeric(groups) is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y=head(y, 10), X=head(X, 10), Z=head(Z, 10), varsel=TRUE, groups=c(1,NA,1,1)), 
               'anyNA(groups) == FALSE is not TRUE', fixed=TRUE)
  ##ztest: only with varsel, numeric, same or fewer columns than Z, no missing
  expect_message(testfit<-kmbayes(y=y, X=X, Z=Z, ztest=c(1,2)), 
                 'ztest should only be specified if varsel=TRUE, resetting varsel to TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y=y, X=X, Z=Z, varsel=TRUE, ztest=c("hello",1)), 
               'is.numeric(ztest) is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y=y, X=X, Z=Z, varsel=TRUE, ztest =c(1,1,1,1,1)), 
               'length(ztest) <= ncol(Z) is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y=y, X=X, Z=Z, varsel=TRUE, ztest =c(1,5,9)), 
               'max(ztest) <= ncol(Z) is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y=y, X=X, Z=Z, varsel=TRUE, ztest=c(1,NA,2)), 
               'anyNA(ztest) == FALSE is not TRUE', fixed=TRUE)
})