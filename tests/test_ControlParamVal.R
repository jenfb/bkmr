
library(bkmr)
library(testthat)
set.seed(111)
dat <- SimData(n = 50, M = 4)
y <- dat$y
Z <- dat$Z
X <- dat$X
test_that("gaussian family", {
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, family="gaussian", control.params=list(a.sigsq=0)), 
               'control.params$a.sigsq > 0 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, family="gaussian", control.params=list(b.sigsq=0)), 
               'control.params$b.sigsq > 0 is not TRUE', fixed=TRUE)
})

test_that("variable selection or not", {
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, varsel=TRUE, control.params=list(a.p0=0)), 
               'control.params$a.p0 > 0 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, varsel=TRUE, control.params=list(b.p0=0)), 
               'control.params$b.p0 > 0 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, varsel=TRUE, control.params=list(r.jump1=0)), 
               'control.params$r.jump1 > 0 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, varsel=TRUE, control.params=list(r.jump2=0)), 
               'control.params$r.jump2 > 0 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, varsel=TRUE, control.params=list(r.muprop=0)), 
               'control.params$r.muprop > 0 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, varsel=FALSE, control.params=list(r.jump=0)), 
               'control.params$r.jump > 0 is not TRUE', fixed=TRUE)
})

myid<-c(seq(1, 48), c(1,2))
test_that("lambda params with and without id", {
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, id=myid, control.params=list(mu.lambda=(1))), 
                 'length(control.params$mu.lambda) == 2 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, id=myid, control.params=list(mu.lambda=c(1,2), sigma.lambda=5)), 
               'length(control.params$sigma.lambda) == 2 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, id=myid, control.params=list(mu.lambda=c(1,2), sigma.lambda=c(5,10), lambda.jump=2)), 
               'length(control.params$lambda.jump) == 2 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, id=myid, control.params=list(mu.lambda=c(1,-1), sigma.lambda=c(5,10), lambda.jump=c(2, 0.2))), 
               'control.params$mu.lambda > 0 are not all TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, id=myid, control.params=list(mu.lambda=c(1,2), sigma.lambda=c(5,-5), lambda.jump=c(2, 0.2))), 
               'control.params$sigma.lambda > 0 are not all TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, id=myid, control.params=list(mu.lambda=c(1,2), sigma.lambda=c(5,10), lambda.jump=c(2, -0.2))), 
               'control.params$lambda.jump > 0 are not all TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, control.params=list(mu.lambda=(-1), sigma.lambda=10, lambda.jump=2)), 
               'control.params$mu.lambda > 0 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, control.params=list(mu.lambda=1, sigma.lambda=(-5), lambda.jump=2)), 
               'control.params$sigma.lambda > 0 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, control.params=list(mu.lambda=1, sigma.lambda=5, lambda.jump=(-2))), 
               'control.params$lambda.jump > 0 is not TRUE', fixed=TRUE)
})
  
test_that("prior distributions", {
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, control.params=list(r.prior="hello")), 
              'rprior == "gamma" | rprior == "unif" | rprior == "invunif" is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, control.params=list(r.prior="gamma", mu.r=0)), 
                 'control.params$mu.r > 0 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, control.params=list(r.prior="gamma", sigma.r=0)), 
                 'control.params$sigma.r > 0 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, control.params=list(r.prior="unif", r.a=(-1))), 
                 'control.params$r.a >= 0 is not TRUE', fixed=TRUE)
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, control.params=list(r.prior="unif", r.b=(-1))), 
                 'control.params$r.b > control.params$r.a is not TRUE', fixed=TRUE)
})
