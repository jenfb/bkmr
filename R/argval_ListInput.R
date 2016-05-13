# Validate control params list
##components of list
##lambda.jump             / default=10
##mu.lambda, sigma.lambda / default=10
##a.p0, b.p0              / default=1
##r.prior                 / default = "gamma", alt=invunif, unif
##a.sigsq, b.sigsq        / default=0.001
##mu.r, sigma.r           / default=5
##r.muprop                / default=1
##r.jump                  / default=0.2
##r.jump1, r.jump2        / default=2, 0.2
##r.a, r.b                / default=0, 100

# 
validateStartingValues <- function(varsel, y, X, Z, starting.values) {
  Ylength <- length(y)
  Xwidth <- ncol(X)
  Zwidth <- ncol(Z) 
  message ("Validating starting.values...")
  ##print(starting.values)
  stopifnot(starting.values$sigsq.eps > 0, starting.values$lambda > 0)
  ##messages only for scalar where vector required, expansion happens in main function
  if (length(starting.values$beta) != Xwidth) {
    message("beta should be a vector of length equal to the number of columns of X.  Input will be repeated or truncated as necessary.")
  }
  ##h.hat length and values
  if (length(starting.values$h.hat) != Ylength) {
    message("h.hat should be a vector of length equal to number of rows in Y.  Input will be repeated or truncated as necessary.")
  }
  for (i in 1:length(starting.values$h.hat)) {
    stopifnot(starting.values$h.hat > 0) 
  }
  ##delta length, 0 or 1 are valid values
  if (length(starting.values$delta) != Zwidth) {
    message("delta should be a vector of length equal to the number of columns of Z.  Input will be repeated or truncated as necessary.")
  }
  for (i in 1:length(starting.values$delta)) {
    stopifnot(starting.values$delta == 1 || starting.values$delta == 0) 
  }
  ## r length depends on varsel, truncate here but expand if necessary in main function
  if (varsel == TRUE) {
    if (length(starting.values$r) != Zwidth) {
      message("r should be a vector of length equal to the number of columns of Z.  Input will be repeated or truncated as necessary.")
    }
  }
  else {
    if (length(starting.values$r) > 1) {
      message("r should a scalar.  Vector input will be truncated.")
      starting.values$r=starting.values$r[1]
    }
  }

  for (i in 1:length(starting.values$r)) {
    stopifnot(starting.values$r > 0) 
  }
}
  
    
validateControlParams <- function(varsel, family, id, control.params) {
  message ("Validating control.params...")
  ##print(control.params)
  if (family == "gaussian"){
    stopifnot(control.params$a.sigsq > 0, control.params$b.sigsq > 0)
  } 
  if (varsel == TRUE) {
    stopifnot(control.params$a.p0 > 0, control.params$b.p0 > 0, control.params$r.jump1 > 0, control.params$r.jump2 > 0, control.params$r.muprop > 0)
  }
  else  {
    stopifnot(control.params$r.jump > 0)
  } ##end varsel-specific stuff
  ##if id, need two elements in mu.lambda, sigma.lambda and lambda.jump
  if (!missing(id)) {
    stopifnot(length(control.params$mu.lambda) == 2, length(control.params$sigma.lambda) == 2, length(control.params$lambda.jump) == 2)
  }
  ##regardless of id, validate each element of these params
  for (i in 1:length(control.params$mu.lambda)) {
    stopifnot(control.params$mu.lambda > 0) 
  }
  for (i in 1:length(control.params$sigma.lambda)) {
    stopifnot(control.params$sigma.lambda > 0) 
  }
  for (i in 1:length(control.params$lambda.jump)) {
    stopifnot(control.params$lambda.jump > 0) 
  }
  rprior=control.params$r.prior
  stopifnot(rprior == "gamma" | rprior == "unif" | rprior == "invunif")
  ##stopifnot(length(intersect (control.params$r.prior, c("gamma","unif","invunif")))>0)
  if (control.params$r.prior == "gamma") {
    stopifnot(control.params$mu.r > 0, control.params$sigma.r > 0)
  }
  else {  
    stopifnot(control.params$r.a >= 0, control.params$r.b > control.params$r.a)
  }
}