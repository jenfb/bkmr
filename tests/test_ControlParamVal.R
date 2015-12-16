source("R/argval_controlParams.R")
validateControlParams(varsel=TRUE, family="gaussian", id, list())
control.params <- list()


library(bkmr)
library(testthat)
set.seed(111)
dat <- SimData(n = 50, M = 4)
y <- dat$y
Z <- dat$Z
X <- dat$X
cp1<-list()
cp1$a.sigsq<-(-1)
cp2<-list()
cp2$b.sigsq<-(-1)
test_that("gaussian family", {
  expect_message(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, family="gaussian", control.params=cp1), 
               'invalid value for a.sigsq, resetting to default .001', fixed=TRUE)
  expect_message(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, family="gaussian", control.params=cp2), 
               'invalid value for b.sigsq, resetting to default .001', fixed=TRUE)
})

cp1<-list()
cp1$a.p0<-(-1)
cp2<-list()
cp2$b.p0<-(-1)
cp3<-list()
cp3$r.jump1<-(-1)
cp4<-list()
cp4$r.jump2<-(-1)
cp5<-list()
cp5$r.muprop<-(-1)
cp6<-list()
cp6$r.jump<-(-1)
test_that("variable selection or not", {
  expect_message(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, varsel=TRUE, control.params=cp1), 
               'invalid value for a.p0 , resetting to default 1', fixed=TRUE)
  expect_message(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, varsel=TRUE, control.params=cp2), 
               'invalid value for b.p0 , resetting to default 1', fixed=TRUE)
  expect_message(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, varsel=TRUE, control.params=cp3), 
               'invalid value for r.jump1, resetting to default 2', fixed=TRUE)
  expect_message(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, varsel=TRUE, control.params=cp4), 
               'invalid value for r.jump2, resetting to default 0.2', fixed=TRUE)
  expect_message(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, varsel=TRUE, control.params=cp5), 
               'invalid value for r.muprop, resetting to default 1', fixed=TRUE)
  expect_message(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, varsel=FALSE, control.params=cp6), 
                 'invalid value for r.jump, resetting to default 0.2', fixed=TRUE)
})

cp1<-list()
cp1$mu.lambda<-(1)
myid<-c(seq(1, 48), c(1,2))
test_that("id", {
  expect_error(testfit<-kmbayes(y = y, Z = Z, X = X, iter = 100, id=myid, control.params=cp6), 
                 'length(control.params$mu.lambda) == 2 & length(control.params$sigma.lambda) ==  .... is not TRUE', fixed=TRUE)
})