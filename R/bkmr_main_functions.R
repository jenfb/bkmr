# makeKpart <- function(r, Z) {
# Kpart <- as.matrix(dist(sqrt(matrix(r, byrow=TRUE, nrow(Z), ncol(Z)))*Z))^2
# Kpart
# }
makeKpart <- function(r, Z1, Z2 = NULL) {
  Z1r <- sweep(Z1, 2, sqrt(r), "*")
  if (is.null(Z2)) {
    Z2r <- Z1r
  } else {
    Z2r <- sweep(Z2, 2, sqrt(r), "*")
  }
  Kpart <- fields::rdist(Z1r, Z2r)^2
  Kpart
}
makeVcomps <- function(r, lambda, Z, data.comps) {
  if (is.null(data.comps$knots)) {
    Kpart <- makeKpart(r, Z)
    V <- diag(1, nrow(Z), nrow(Z)) + lambda[1]*exp(-Kpart)
    if (data.comps$nlambda == 2) {
      V <- V + lambda[2]*data.comps$crossTT
    }
    cholV <- chol(V)
    Vinv <- chol2inv(cholV)
    logdetVinv <- -2*sum(log(diag(cholV)))
    Vcomps <- list(Vinv = Vinv, logdetVinv = logdetVinv)
  } else {## predictive process approach
    ## note: currently does not work with random intercept model
    nugget <- 0.001
    n0 <- nrow(Z)
    n1 <- nrow(data.comps$knots)
    nall <- n0 + n1
    # Kpartall <- makeKpart(r, rbind(Z, data.comps$knots))
    # Kall <- exp(-Kpartall)
    # K0 <- Kall[1:n0, 1:n0 ,drop=FALSE]
    # K1 <- Kall[(n0+1):nall, (n0+1):nall ,drop=FALSE]
    # K10 <- Kall[(n0+1):nall, 1:n0 ,drop=FALSE]
    K1 <- exp(-makeKpart(r, data.comps$knots))
    K10 <- exp(-makeKpart(r, data.comps$knots, Z))
    Q <- K1 + diag(nugget, n1, n1)
    R <- Q + lambda[1]*tcrossprod(K10)
    cholQ <- chol(Q)
    cholR <- chol(R)
    Qinv <- chol2inv(cholQ)
    Rinv <- chol2inv(cholR)
    Vinv <- diag(1, n0, n0) - lambda[1]*t(K10) %*% Rinv %*% K10
    logdetVinv <- 2*sum(log(diag(cholQ))) - 2*sum(log(diag(cholR)))
    Vcomps <- list(Vinv = Vinv, logdetVinv = logdetVinv, cholR = cholR, Q = Q, K10 = K10, Qinv = Qinv, Rinv = Rinv)
  }
  Vcomps
}

#' Fit Bayesian kernel machine regression
#'
#' Fits the Bayesian kernel machine regression (BKMR) model using Markov chain Monte Carlo (MCMC) methods.
#'
#' @export
#'
#' @param y a vector of outcome data of length \code{n}.
#' @param Z an \code{n}-by-\code{M} matrix of predictor variables to be included in the \code{h} function. Each row represents an observation and each column represents an predictor.
#' @param X an \code{n}-by-\code{K} matrix of covariate data where each row represents an observation and each column represents a covariate. Should not contain an intercept column.
#' @param iter number of iterations to run the sampler
#' @param family a description of the error distribution and link function to be used in the model. Currently implemented for \code{gaussian} and \code{binomial} families.
#' @param id optional vector (of length \code{n}) of grouping factors for fitting a model with a random intercept. If NULL then no random intercept will be included.
#' @param verbose TRUE or FALSE: flag indicating whether to print intermediate diagnostic information during the model fitting.
#' @param Znew optional matrix of new predictor values at which to predict \code{h}, where each row represents a new observation. This will slow down the model fitting, and can be done as a post-processing step using \code{\link{SamplePred}}
#' @param starting.values list of starting values for each parameter. If not specified default values will be chosen.
#' @param control.params list of parameters specifying the prior distributions and tuning parameters for the MCMC algorithm. If not specified default values will be chosen.
#' @param varsel TRUE or FALSE: indicator for whether to conduct variable selection on the Z variables in \code{h}
#' @param groups optional vector (of length \code{M}) of group indictors for fitting hierarchical variable selection if varsel=TRUE. If varsel=TRUE without group specification, component-wise variable selections will be performed.
#' @param knots optional matrix of knot locations for implementing the Gaussian predictive process of Banerjee et al (2008). Currently only implemented for models without a random intercept.
#' @param ztest optional vector indicating on which variables in Z to conduct variable selection (the remaining variables will be forced into the model).
#' @param rmethod for those predictors being forced into the \code{h} function, the method for sampling the \code{r[m]} values. Takes the value of 'varying' to allow separate \code{r[m]} for each predictor; 'equal' to force the same \code{r[m]} for each predictor; or 'fixed' to fix the \code{r[m]} to their starting values
#' @param est.h TRUE or FALSE: indicator for whether to sample from the posterior distribution of the subject-specific effects h_i within the main sampler. This will slow down the model fitting.
#' @return an object of class "bkmrfit", which has the associated methods:
#' \itemize{
#'   \item \code{\link{print}} (i.e., \code{\link{print.bkmrfit}}) 
#'   \item \code{\link{summary}} (i.e., \code{\link{summary.bkmrfit}})
#' }
#' 
#' @seealso For guided examples, go to \url{https://jenfb.github.io/bkmr/overview.html}
#' @references Bobb, JF, Valeri L, Claus Henn B, Christiani DC, Wright RO, Mazumdar M, Godleski JJ, Coull BA (2015). Bayesian Kernel Machine Regression for Estimating the Health Effects of Multi-Pollutant Mixtures. Biostatistics 16, no. 3: 493-508.
#' @references Banerjee S, Gelfand AE, Finley AO, Sang H (2008). Gaussian predictive process models for large spatial data sets. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70(4), 825-848.
#' @import utils
kmbayes <- function(y, Z, X = NULL, iter = 1000, family = "gaussian", id = NULL, verbose = TRUE, Znew = NULL, starting.values = NULL, control.params = NULL, varsel = FALSE, groups = NULL, knots = NULL, ztest = NULL, rmethod = "varying", est.h = FALSE) {
  
  missingX <- is.null(X)
  if (missingX) X <- matrix(0, length(y), 1)
  hier_varsel <- !is.null(groups)
  
  ##Argument check 1, required arguments without defaults
  ##check vector/matrix sizes
  stopifnot (length(y) > 0, is.numeric(y), anyNA(y) == FALSE)
  if (inherits(class(Z), "matrix") == FALSE)  Z <- as.matrix(Z)
  stopifnot (is.numeric(Z), nrow(Z) == length(y), anyNA(Z) == FALSE)
  if (inherits(class(X), "matrix") == FALSE)  X <- as.matrix(X)
  stopifnot (is.numeric(X), nrow(X) == length(y), anyNA(X) == FALSE) 
  
  ##Argument check 2: for those with defaults, write message and reset to default if invalid
  if (iter < 1) {
    message ("invalid input for iter, resetting to default value 1000")
    nsamp <- 1000
  } else {
    nsamp <- iter
  }
  if (!family %in% c("gaussian", "binomial")) {
    stop("family", family, "not yet implemented; must specify either 'gaussian' or 'binomial'")
  }
  if (family == "binomial") {
    message("Fitting probit regression model")
    if (!all(y %in% c(0, 1))) {
      stop("When family == 'binomial', y must be a vector containing only zeros and ones")
    }
  }
  if (rmethod != "varying" & rmethod != "equal" & rmethod != "fixed") {
    message ("invalid value for rmethod, resetting to default varying")
    rmethod <- "varying"
  }
  if (verbose != FALSE & verbose != TRUE) {
    message ("invalid value for verbose, resetting to default FALSE")
    verbose <- FALSE
  }
  if (varsel != FALSE & varsel != TRUE) {
    message ("invalid value for varsel, resetting to default FALSE")
    varsel <- FALSE
  }
  
  ##Argument check 3: the rest id (below) znew, knots, groups, ztest
  if (!is.null(id)) { 
    stopifnot(length(id) == length(y), anyNA(id) == FALSE)
    if (!is.null(knots)) { 
      message ("knots cannot be specified with id, resetting knots to null")
      knots<-NA
    }
  }
  if (!is.null(Znew)) { 
    if (class(Znew) != "matrix")  Znew <- as.matrix(Znew)
    stopifnot(is.numeric(Znew), ncol(Znew) == ncol(Z), anyNA(Znew) == FALSE)
  }
  if (!is.null(knots)) { 
    if (class(knots) != "matrix")  knots <- as.matrix(knots)
    stopifnot(is.numeric(knots), ncol(knots )== ncol(Z), anyNA(knots) == FALSE)
  }
  if (!is.null(groups)) { 
    if (varsel == FALSE) {
      message ("groups should only be specified if varsel=TRUE, resetting varsel to TRUE")
      varsel <- TRUE
    } else {
      stopifnot(is.numeric(groups), length(groups) == ncol(Z), anyNA(groups) == FALSE)
    }
  }
  if (!is.null(ztest)) { 
    if (varsel == FALSE) {
      message ("ztest should only be specified if varsel=TRUE, resetting varsel to TRUE")
      varsel <- TRUE
    } else {
      stopifnot(is.numeric(ztest), length(ztest) <= ncol(Z), anyNA(ztest) == FALSE, max(ztest) <= ncol(Z) )
    }
  }
  
  ## start JB code
  if (!is.null(id)) { ## for random intercept model
    randint <- TRUE
    id <- as.numeric(as.factor(id))
    nid <- length(unique(id))
    nlambda <- 2
    
    ## matrix that multiplies the random intercept vector
    TT <- matrix(0, length(id), nid)
    for (i in 1:nid) {
      TT[which(id == i), i] <- 1
    }
    crossTT <- tcrossprod(TT)
    rm(TT, nid)
  } else {
    randint <- FALSE
    nlambda <- 1
    crossTT <- 0
  }
  data.comps <- list(randint = randint, nlambda = nlambda, crossTT = crossTT)
  if (!is.null(knots)) data.comps$knots <- knots
  rm(randint, nlambda, crossTT)
  
  ## create empty matrices to store the posterior draws in
  chain <- list(h.hat = matrix(0, nsamp, nrow(Z)),
                beta = matrix(0, nsamp, ncol(X)),
                lambda = matrix(NA, nsamp, data.comps$nlambda),
                sigsq.eps = rep(NA, nsamp),
                r = matrix(NA, nsamp, ncol(Z)),
                acc.r = matrix(0, nsamp, ncol(Z)),
                acc.lambda = matrix(0, nsamp, data.comps$nlambda),
                delta = matrix(1, nsamp, ncol(Z))
  )
  if (varsel) {
    chain$acc.rdelta <- rep(0, nsamp)
    chain$move.type <- rep(0, nsamp)
  }
  if (family == "binomial") {
    chain$ystar <- matrix(0, nsamp, length(y))
  }
  
  ## components to predict h(Znew)
  if (!is.null(Znew)) {
    if (is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
    if (class(Znew) == "data.frame") Znew <- data.matrix(Znew)
    if (ncol(Z) != ncol(Znew)) {
      stop("Znew must have the same number of columns as Z")
    }
    ##Kpartall <- as.matrix(dist(rbind(Z,Znew)))^2
    chain$hnew <- matrix(0,nsamp,nrow(Znew))
    colnames(chain$hnew) <- rownames(Znew)
  }
  
  ## components if model selection is being done
  if (varsel) {
    if (is.null(ztest)) {
      ztest <- 1:ncol(Z)
    }
    rdelta.update <- rdelta.comp.update
  } else {
    ztest <- NULL
  }
  
  ## control parameters
  control.params.default <- list(lambda.jump = rep(10, data.comps$nlambda), mu.lambda = rep(10, data.comps$nlambda), sigma.lambda = rep(10, data.comps$nlambda), a.p0 = 1, b.p0 = 1, r.prior = "invunif", a.sigsq = 1e-3, b.sigsq = 1e-3, mu.r = 5, sigma.r = 5, r.muprop = 1, r.jump = 0.1, r.jump1 = 2, r.jump2 = 0.1, r.a = 0, r.b = 100)
  if (!is.null(control.params)){
    control.params <- modifyList(control.params.default, as.list(control.params))
    validateControlParams(varsel, family, id, control.params)
  } else {
    control.params <- control.params.default
  }
  
  control.params$r.params <- with(control.params, list(mu.r = mu.r, sigma.r = sigma.r, r.muprop = r.muprop, r.jump = r.jump, r.jump1 = r.jump1, r.jump2 = r.jump2, r.a = r.a, r.b = r.b))
  
  ## components if grouped model selection is being done
  if (!is.null(groups)) {
    if (!varsel) {
      stop("if doing grouped variable selection, must set varsel = TRUE")
    }
    rdelta.update <- rdelta.group.update
    control.params$group.params <- list(groups = groups, sel.groups = sapply(unique(groups), function(x) min(seq_along(groups)[groups == x])), neach.group = sapply(unique(groups), function(x) sum(groups %in% x)))
  }
  
  ## specify functions for doing the Metropolis-Hastings steps to update r
  e <- environment()
  rfn <- set.r.MH.functions(r.prior = control.params$r.prior)
  rprior.logdens <- rfn$rprior.logdens
  environment(rprior.logdens) <- e
  rprop.gen1 <- rfn$rprop.gen1
  environment(rprop.gen1) <- e
  rprop.logdens1 <- rfn$rprop.logdens1
  environment(rprop.logdens1) <- e
  rprop.gen2 <- rfn$rprop.gen2
  environment(rprop.gen2) <- e
  rprop.logdens2 <- rfn$rprop.logdens2
  environment(rprop.logdens2) <- e
  rprop.gen <- rfn$rprop.gen
  environment(rprop.gen) <- e
  rprop.logdens <- rfn$rprop.logdens
  environment(rprop.logdens) <- e
  rm(e, rfn)
  
  ## initial values
  starting.values0 <- list(h.hat = 1, beta = NULL, sigsq.eps = NULL, r = 1, lambda = 10, delta = 1)
  if (is.null(starting.values)) {
    starting.values <- starting.values0
  } else {
    starting.values <- modifyList(starting.values0, starting.values)
    validateStartingValues (varsel, y, X, Z, starting.values)
  }
  if (family == "gaussian") {
    if (is.null(starting.values$beta) | is.null(starting.values$sigsq.eps)) {
      lmfit0 <- lm(y ~ Z + X)
      if (is.null(starting.values$beta)) {
        coefX <- coef(lmfit0)[grep("X", names(coef(lmfit0)))]
        starting.values$beta <- unname(ifelse(is.na(coefX), 0, coefX))
      }
      if (is.null(starting.values$sigsq.eps)) {
        starting.values$sigsq.eps <- summary(lmfit0)$sigma^2
      }
    } 
  } else if (family == "binomial") {
    starting.values$sigsq.eps <- 1 ## always equal to 1
    if (is.null(starting.values$beta) | is.null(starting.values$ystar)) {
      probitfit0 <- try(glm(y ~ Z + X, family = binomial(link = "probit")))
      if (!inherits(probitfit0, "try-error")) {
        if (is.null(starting.values$beta)) {
          coefX <- coef(probitfit0)[grep("X", names(coef(probitfit0)))]
          starting.values$beta <- unname(ifelse(is.na(coefX), 0, coefX))
        }
        if (is.null(starting.values$ystar)) {
          #prd <- predict(probitfit0)
          #starting.values$ystar <- ifelse(y == 1, abs(prd), -abs(prd))
          starting.values$ystar <- ifelse(y == 1, 1/2, -1/2)
        }
      } else {
        starting.values$beta <- 0
        starting.values$ystar <- ifelse(y == 1, 1/2, -1/2)
      }
    } 
  }
    
  ##print (starting.values)
  ##truncate vectors that are too long
  if (length(starting.values$h.hat) > length(y)) {
    starting.values$h.hat <- starting.values$h.hat[1:length(y)]
  }
  if (length(starting.values$beta) > ncol(X)) {
    starting.values$beta <- starting.values$beta[1:ncol(X)]
  }
  if (length(starting.values$delta) > ncol(Z)) {
    starting.values$delta <- starting.values$delta[1:ncol(Z)]
  }
  if (varsel==FALSE & length(starting.values$r) > 1) {
    starting.values$r <- starting.values$r[1]
  } else if (length(starting.values$r) > ncol(Z)) {
    starting.values$r <- starting.values$r[1:ncol(Z)]
  }

  chain$h.hat[1, ] <- starting.values$h.hat
  chain$beta[1, ] <- starting.values$beta
  chain$lambda[1, ] <- starting.values$lambda
  chain$sigsq.eps[1] <- starting.values$sigsq.eps
  chain$r[1, ] <- starting.values$r
  if (varsel) {
    chain$delta[1,ztest] <- starting.values$delta
  }
  if (family == "binomial") {
    chain$ystar[1, ] <- starting.values$ystar
    chain$sigsq.eps[] <- starting.values$sigsq.eps ## does not get updated
  }
  if (!is.null(groups)) {
    ## make sure starting values are consistent with structure of model
    if (!all(sapply(unique(groups), function(x) sum(chain$delta[1, ztest][groups == x])) == 1)) {
      # warning("Specified starting values for delta not consistent with model; using default")
      starting.values$delta <- rep(0, length(groups))
      starting.values$delta[sapply(unique(groups), function(x) min(which(groups == x)))] <- 1
    }
    chain$delta[1,ztest] <- starting.values$delta
    chain$r[1,ztest] <- ifelse(chain$delta[1,ztest] == 1, chain$r[1,ztest], 0)
  }
  chain$est.h <- est.h
  
  ## components
  Vcomps <- makeVcomps(r = chain$r[1, ], lambda = chain$lambda[1, ], Z = Z, data.comps = data.comps)

  ## start sampling ####
  chain$time1 <- Sys.time()
  for (s in 2:nsamp) {

    ## continuous version of outcome (latent outcome under binomial probit model)
    if (family == "gaussian") {
      ycont <- y
    } else if (family == "binomial") {
      if (est.h) {
        chain$ystar[s,] <- ystar.update(y = y, X = X, beta = chain$beta[s - 1,], h = chain$h[s - 1, ])
      } else {
        chain$ystar[s,] <- ystar.update.noh(y = y, X = X, beta = chain$beta[s - 1,], Vinv = Vcomps$Vinv, ystar = chain$ystar[s - 1, ])
      }
      ycont <- chain$ystar[s, ]
    }
    
    ## generate posterior samples from marginalized distribution P(beta, sigsq.eps, lambda, r | y)
    
    ## beta
    if (!missingX) {
      chain$beta[s,] <- beta.update(X = X, Vinv = Vcomps$Vinv, y = ycont, sigsq.eps = chain$sigsq.eps[s - 1])
    }
      
    ## \sigma_\epsilon^2
    if (family == "gaussian") {
      chain$sigsq.eps[s] <- sigsq.eps.update(y = ycont, X = X, beta = chain$beta[s,], Vinv = Vcomps$Vinv, a.eps = control.params$a.sigsq, b.eps = control.params$b.sigsq)
    }
    
    ## lambda
    lambdaSim <- chain$lambda[s - 1,]
    for (comp in 1:data.comps$nlambda) {
      varcomps <- lambda.update(r = chain$r[s - 1,], delta = chain$delta[s - 1,], lambda = lambdaSim, whichcomp = comp, y = ycont, X = X, Z = Z, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, data.comps = data.comps, control.params = control.params)
      lambdaSim <- varcomps$lambda
      if (varcomps$acc) {
        Vcomps <- varcomps$Vcomps
        chain$acc.lambda[s,comp] <- varcomps$acc
      }
    }
    chain$lambda[s,] <- lambdaSim
    
    ## r
    rSim <- chain$r[s - 1,]
    comp <- which(!1:ncol(Z) %in% ztest)
    if (length(comp) != 0) {
      if (rmethod == "equal") { ## common r for those variables not being selected
        varcomps <- r.update(r = rSim, whichcomp = comp, delta = chain$delta[s - 1,], lambda = chain$lambda[s,], y = ycont, X = X, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z, data.comps = data.comps, control.params = control.params, rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1, rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2, rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen, rprop.logdens = rprop.logdens)
        rSim <- varcomps$r
        if (varcomps$acc) {
          Vcomps <- varcomps$Vcomps
          chain$acc.r[s, comp] <- varcomps$acc
        }
      } else if (rmethod == "varying") { ## allow a different r_m
        for (whichr in comp) {
          varcomps <- r.update(r = rSim, whichcomp = whichr, delta = chain$delta[s - 1,], lambda = chain$lambda[s,], y = ycont, X = X, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z, data.comps = data.comps, control.params = control.params, rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1, rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2, rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen, rprop.logdens = rprop.logdens)
          rSim <- varcomps$r
          if (varcomps$acc) {
            Vcomps <- varcomps$Vcomps
            chain$acc.r[s, whichr] <- varcomps$acc
          }
        }
      }
    }
    ## for those variables being selected: joint posterior of (r,delta)
    if (varsel) {
      varcomps <- rdelta.update(r = rSim, delta = chain$delta[s - 1,], lambda = chain$lambda[s,], y = ycont, X = X, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z, ztest = ztest, data.comps = data.comps, control.params = control.params, rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1, rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2, rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen, rprop.logdens = rprop.logdens)
      chain$delta[s,] <- varcomps$delta
      rSim <- varcomps$r
      chain$move.type[s] <- varcomps$move.type
      if (varcomps$acc) {
        Vcomps <- varcomps$Vcomps
        chain$acc.rdelta[s] <- varcomps$acc
      }
    }
    chain$r[s,] <- rSim
    
    ###################################################
    ## generate posterior sample of h(z) from its posterior P(h | beta, sigsq.eps, lambda, r, y)
    
    if (est.h) {
      hcomps <- h.update(lambda = chain$lambda[s,], Vcomps = Vcomps, sigsq.eps = chain$sigsq.eps[s], y = ycont, X = X, beta = chain$beta[s,], r = chain$r[s,], Z = Z, data.comps = data.comps)
      chain$h.hat[s,] <- hcomps$hsamp
      if (!is.null(hcomps$hsamp.star)) { ## GPP
        Vcomps$hsamp.star <- hcomps$hsamp.star
      }
      rm(hcomps)
    }
      
    ###################################################
    ## generate posterior samples of h(Znew) from its posterior P(hnew | beta, sigsq.eps, lambda, r, y)
    
    if (!is.null(Znew)) {
      chain$hnew[s,] <- newh.update(Z = Z, Znew = Znew, Vcomps = Vcomps, lambda = chain$lambda[s,], sigsq.eps = chain$sigsq.eps[s], r = chain$r[s,], y = ycont, X = X, beta = chain$beta[s,], data.comps = data.comps)
    }
    
    ###################################################
    ## print details of the model fit so far
    opts <- set_verbose_opts(
      verbose_freq = control.params$verbose_freq, 
      verbose_digits = control.params$verbose_digits,
      verbose_show_ests = control.params$verbose_show_ests
      )
    print_diagnostics(verbose = verbose, opts = opts, curr_iter = s, tot_iter = nsamp, chain = chain, varsel = varsel, hier_varsel = hier_varsel, ztest = ztest, Z = Z, groups = groups)
   
  }
  control.params$r.params <- NULL
  chain$time2 <- Sys.time()
  chain$iter <- nsamp
  chain$family <- family
  chain$starting.values <- starting.values
  chain$control.params <- control.params
  chain$X <- X
  chain$Z <- Z
  chain$y <- y
  chain$ztest <- ztest
  chain$data.comps <- data.comps
  if (!is.null(Znew)) chain$Znew <- Znew
  if (!is.null(groups)) chain$groups <- groups
  chain$varsel <- varsel
  class(chain) <- c("bkmrfit", class(chain))
  chain
}

#' Print basic summary of BKMR model fit
#'
#' \code{print} method for class "bkmrfit"
#'
#' @param x an object of class "bkmrfit"
#' @param digits the number of digits to show when printing
#' @param ...	further arguments passed to or from other methods.
#'  
#' @export
print.bkmrfit <- function(x, digits = 5, ...) {
  cat("Fitted object of class 'bkmrfit'\n")
  cat("Iterations:", x$iter, "\n")
  cat("Outcome family:", x$family, ifelse(x$family == "binomial", "(probit link)", ""), "\n")
  cat("Model fit on:", as.character(x$time2), "\n")
}

#' Summarizing BKMR model fits
#'
#' \code{summary} method for class "bkmrfit"
#'
#' @param object an object of class "bkmrfit"
#' @param q quantiles of posterior distribution to show
#' @param digits the number of digits to show when printing
#' @param show_ests logical; if \code{TRUE}, prints summary statistics of posterior distribution
#' @param show_MH logical; if \code{TRUE}, prints acceptance rates from the Metropolis-Hastings algorithm
#' @param ...	further arguments passed to or from other methods.
#'  
#' @export
summary.bkmrfit <- function(object, q = c(0.025, 0.975), digits = 5, show_ests = TRUE, show_MH = TRUE, ...) {
  x <- object
  elapsed_time <- difftime(x$time2, x$time1)

  print(x, digits = digits)  
  cat("Running time: ", round(elapsed_time, digits), attr(elapsed_time, "units"), "\n")
  cat("\n")
  
  if (show_MH) {
    cat("Acceptance rates for Metropolis-Hastings algorithm:\n")
    accep_rates <- data.frame()
    ## lambda
    nm <- "lambda"
    rate <- colMeans(x$acc.lambda[2:x$iter, ,drop = FALSE])
    if (length(rate) > 1) nm <- paste0(nm, seq_along(rate))
    accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
    ## r_m
    if (!x$varsel) {
      nm <- "r"
      rate <- colMeans(x$acc.r[2:x$iter, , drop = FALSE])
      if (length(rate) > 1) nm <- paste0(nm, seq_along(rate))
      accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
    } else {
      nm <- "r/delta (overall)"
      rate <- mean(x$acc.rdelta[2:x$iter])
      accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      ##
      nm <- "r/delta  (move 1)"
      rate <- mean(x$acc.rdelta[2:x$iter][x$move.type[2:x$iter] == 1])
      accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      ##
      nm <- "r/delta  (move 2)"
      rate <- mean(x$acc.rdelta[2:x$iter][x$move.type[2:x$iter] == 2])
      accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      if (!is.null(x$groups)) {
        nm <- "r/delta  (move 3)"
        rate <- mean(x$acc.rdelta[2:x$iter][x$move.type[2:x$iter] == 3])
        accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      }
    }
    print(accep_rates)
  }
  if (show_ests) {
    sel <- with(x, seq(floor(iter/2) + 1, iter))
    cat("\nParameter estimates (based on iterations ", min(sel), "-", max(sel), "):\n", sep = "")
    ests <- ExtractEsts(x, q = q, sel = sel)
    if (!is.null(ests$h)) {
      ests$h <- ests$h[c(1,2,nrow(ests$h)), ]
    }
    if (!is.null(ests$ystar)) {
      ests$ystar <- ests$ystar[c(1,2,nrow(ests$ystar)), ]
    }
    summ <- with(ests, rbind(beta, sigsq.eps, r, lambda))
    if (!is.null(ests$h)) {
      summ <- rbind(summ, ests$h)
    }
    if (!is.null(ests$ystar)) {
      summ <- rbind(summ, ests$ystar)
    }
    summ <- data.frame(param = rownames(summ), round(summ, digits))
    rownames(summ) <- NULL
    print(summ)
    if (x$varsel) {
      cat("\nPosterior inclusion probabilities:\n")
      pips <- ExtractPIPs(x)
      pips[, -1] <- round(pips[, -1], digits)
      print(pips)
    }
  }
}
