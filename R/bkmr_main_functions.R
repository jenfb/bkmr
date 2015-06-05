# makeKpart <- function(r, Z) {
	# Kpart <- as.matrix(dist(sqrt(matrix(r, byrow=TRUE, nrow(Z), ncol(Z)))*Z))^2
	# Kpart
# }
makeKpart <- function(r, Z1, Z2) {
	Z1r <- sweep(Z1, 2, sqrt(r), "*")
	if(missing(Z2)) {
		Z2r <- Z1r
	} else {
		Z2r <- sweep(Z2, 2, sqrt(r), "*")
	}
	Kpart <- fields::rdist(Z1r, Z2r)^2
	Kpart
}
makeVcomps <- function(r, lambda, Z, data.comps) {
	if(is.null(data.comps$knots)) {
		Kpart <- makeKpart(r, Z)
		V <- diag(1, nrow(Z), nrow(Z)) + lambda[1]*exp(-Kpart)
		if(data.comps$nlambda == 2) {
			V <- V + lambda[2]*data.comps$crossTT
		}
		cholV <- chol(V)
		Vinv <- chol2inv(cholV)
		logdetVinv <- -2*sum(log(diag(cholV)))
		Vcomps <- list(Vinv = Vinv, logdetVinv = logdetVinv)
	} else { ## predictive process approach
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

#' @export
#' @param equalr Set to \code{TRUE} if you want to allow each exposure variable to have a separate smoothing parameter
kmbayes <- function(y, covar, expos, iter = 1000, id, quiet=TRUE, Znew, starting.values=list(), control.params=list(), modsel=FALSE, groups, knots, ztest, equalr = FALSE) {

    X <- covar
    Z <- expos
    nsamp <- iter

	if(!missing(id)) { ## for random intercept model
		randint <- TRUE
		id <- as.numeric(as.factor(id))
		nid <- length(unique(id))
		nlambda <- 2

		## matrix that multiplies the random intercept vector
		TT <- matrix(0, length(id), nid)
		for(i in 1:nid) {
			TT[which(id==i),i] <- 1
		}
		crossTT <- tcrossprod(TT)
		rm(TT, nid)
	} else {
		randint <- FALSE
		nlambda <- 1
		crossTT <- 0
	}
	data.comps <- list(randint = randint, nlambda = nlambda, crossTT = crossTT)
	if(!missing(knots)) data.comps$knots <- knots
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
	if(modsel) {
		chain$acc.rdelta <- rep(0, nsamp)
		chain$move.type <- rep(0, nsamp)
	}

	## components to predict h(Znew)
	if(!missing(Znew)) {
		if(is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
		if(class(Znew) == "data.frame") Znew <- data.matrix(Znew)
		if(ncol(Z) != ncol(Znew)) {
			stop("Znew must have the same number of columns as Z")
		}
		##Kpartall <- as.matrix(dist(rbind(Z,Znew)))^2
		chain$hnew <- matrix(0,nsamp,nrow(Znew))
		colnames(chain$hnew) <- rownames(Znew)
	}

	## components if model selection is being done
	if(modsel) {
		if(missing(ztest)) {
			ztest <- 1:ncol(Z)
		}
		rdelta.update <- rdelta.comp.update
	} else {
		ztest <- NULL
	}

	## control parameters
	control.params <- modifyList(list(lambda.jump=10, mu.lambda=10, sigma.lambda=10, a.p0=1, b.p0=1, r.prior="gamma", a.sigsq=1e-3, b.sigsq=1e-3, r.params=list(mu.r=5, sigma.r=5, r.muprop=1, r.jump=0.2, r.jump1 = 2, r.jump2 = 0.2)), control.params)

	## components if grouped model selection is being done
	if(!missing(groups)) {
		if(!modsel) {
			stop("if doing grouped variable selection, must set modsel = TRUE")
		}
		rdelta.update <- rdelta.group.update
		control.params$group.params <- list(groups = groups, sel.groups = sapply(unique(groups), function(x) min(seq_along(groups)[groups == x])), neach.group = sapply(unique(groups), function(x) sum(groups %in% x)))
	}

	## specify functions for doing the Metropolis-Hastings steps to update r
    e <- environment()
    rfn <- set.r.MH.functions(r.prior = control.params$r.prior)
	assign("rprior.logdens", rfn$rprior.logdens, envir = parent.env(e))
	assign("rprop.gen1", rfn$rprop.gen1, envir = parent.env(e))
	assign("rprop.logdens1", rfn$rprop.logdens1, envir = parent.env(e))
	assign("rprop.gen2", rfn$rprop.gen2, envir = parent.env(e))
	assign("rprop.logdens2", rfn$rprop.logdens2, envir = parent.env(e))
	assign("rprop.gen", rfn$rprop.gen, envir = parent.env(e))
	assign("rprop.logdens", rfn$rprop.logdens, envir = parent.env(e))
    rm(e, rfn)

	## initial values
    if (is.null(starting.values$beta) | is.null(starting.values$sigsq.eps)) {
        lmfit0 <- lm(y ~ expos + covar)
        if (is.null(starting.values$beta)) {
            starting.values$beta <- coef(lmfit0)[grep("covar", names(coef(lmfit0)))]
        }
        if (is.null(starting.values$sigsq.eps)) {
            starting.values$sigsq.eps <- summary(lmfit0)$sigma^2
        }
    }
	starting.values <- modifyList(list(h.hat=1, beta=1, sigsq.eps=1, r=1, lambda=10, delta=1), starting.values)
	chain$h.hat[1,] <- starting.values$h.hat
	chain$beta[1,] <- starting.values$beta
	chain$lambda[1,] <- starting.values$lambda
	chain$sigsq.eps[1] <- starting.values$sigsq.eps
	chain$r[1,] <- starting.values$r
	if(modsel) {
		chain$delta[1,ztest] <- starting.values$delta
	}
	if(!missing(groups)) {
		## make sure starting values are consistent with structure of model
		if(!all(sapply(unique(groups), function(x) sum(chain$delta[1,ztest][groups == x])) == 1)) {
			# warning("Specified starting values for delta not consistent with model; using default")
			starting.values$delta <- rep(0, length(groups))
			starting.values$delta[sapply(unique(groups), function(x) min(which(groups == x)))] <- 1
		}
		chain$delta[1,ztest] <- starting.values$delta
		chain$r[1,ztest] <- ifelse(chain$delta[1,ztest] == 1, chain$r[1,ztest], 0)
	}

	## components
	Vcomps <- makeVcomps(r = chain$r[1, ], lambda = chain$lambda[1, ], Z = Z, data.comps = data.comps)

	chain$time1 <- Sys.time()
	for(s in 2:nsamp) {
		###################################################
		## generate posterior samples from marginalized distribution P(beta, sigsq.eps, lambda, r | y)

		## beta
		chain$beta[s,] <- beta.update(X = X, Vinv = Vcomps$Vinv, y = y, sigsq.eps = chain$sigsq.eps[s-1])

		## \sigma_\epsilon^2
		chain$sigsq.eps[s] <- sigsq.eps.update(y = y, X = X, beta = chain$beta[s,], Vinv = Vcomps$Vinv, a.eps = control.params$a.sigsq, b.eps = control.params$b.sigsq)

		## lambda
		lambdaSim <- chain$lambda[s-1,]
		for(comp in 1:data.comps$nlambda) {
			varcomps <- lambda.update(r = chain$r[s-1,], delta = chain$delta[s-1,], lambda = lambdaSim, whichcomp = comp, y = y, X = X, Z = Z, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, data.comps = data.comps, control.params = control.params)
			lambdaSim <- varcomps$lambda
			if(varcomps$acc) {
				Vcomps <- varcomps$Vcomps
				chain$acc.lambda[s,comp] <- varcomps$acc
			}
		}
		chain$lambda[s,] <- lambdaSim

		## r
        rSim <- chain$r[s-1,]
        comp <- which(!1:ncol(Z) %in% ztest)
		if(length(comp) != 0) {
            if(equalr) { ## common r for those variables not being selected
                varcomps <- r.update(r = rSim, whichcomp = comp, delta = chain$delta[s-1,], lambda = chain$lambda[s,], y = y, X = X, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z, data.comps = data.comps, control.params = control.params)
                rSim <- varcomps$r
                if(varcomps$acc) {
                    Vcomps <- varcomps$Vcomps
                    chain$acc.r[s, comp] <- varcomps$acc
                }
            } else { ## allow a different r_m
                for(whichr in comp) {
                    varcomps <- r.update(r = rSim, whichcomp = whichr, delta = chain$delta[s-1,], lambda = chain$lambda[s,], y = y, X = X, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z, data.comps = data.comps, control.params = control.params)
                    rSim <- varcomps$r
                    if(varcomps$acc) {
                        Vcomps <- varcomps$Vcomps
                        chain$acc.r[s, whichr] <- varcomps$acc
                    }
                }
            }
		}
		## for those variables being selected: joint posterior of (r,delta)
		if(modsel) {
			varcomps <- rdelta.update(r = rSim, delta = chain$delta[s-1,], lambda = chain$lambda[s,], y = y, X = X, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z, ztest = ztest, data.comps = data.comps, control.params = control.params)
			chain$delta[s,] <- varcomps$delta
			rSim <- varcomps$r
			chain$move.type[s] <- varcomps$move.type
			if(varcomps$acc) {
				Vcomps <- varcomps$Vcomps
				chain$acc.rdelta[s] <- varcomps$acc
			}
		}
		chain$r[s,] <- rSim

		###################################################
		## generate posterior sample of h(z) from its posterior P(h | beta, sigsq.eps, lambda, r, y)

		hcomps <- h.update(lambda = chain$lambda[s,], Vcomps = Vcomps, sigsq.eps = chain$sigsq.eps[s], y = y, X = X, beta = chain$beta[s,], r = chain$r[s,], Z = Z)
		chain$h.hat[s,] <- hcomps$hsamp
		if(!is.null(hcomps$hsamp.star)) { ## GPP
			Vcomps$hsamp.star <- hcomps$hsamp.star
		}
		rm(hcomps)

		###################################################
		## generate posterior samples of h(Znew) from its posterior P(hnew | beta, sigsq.eps, lambda, r, y)

		if(!missing(Znew)) {
			chain$hnew[s,] <- newh.update(Z = Z, Znew = Znew, Vcomps = Vcomps, lambda = chain$lambda[s,], sigsq.eps = chain$sigsq.eps[s], r = chain$r[s,], y = y, X = X, beta = chain$beta[s,], data.comps = data.comps)
		}

		###################################################
		## print details of the model fit so far
		if(s%%(nsamp/10)==0 & !quiet) {
			print(s)
			cat(round(colMeans(chain$acc.lambda[1:s, ,drop=FALSE]),4), "   lam accept rate\n")
			cat(round(colMeans(chain$acc.r[2:s, ]),4), "   r nosel accept rate\n")
			if(modsel) {
				cat(round(mean(chain$acc.rdelta[2:s]),4), "   rdelt accept rate\n")
				cat(round(mean(chain$acc.rdelta[2:s][chain$move.type[2:s] == 1]),4), "   rdelt[move 1] accept rate\n")
				cat(round(mean(chain$acc.rdelta[2:s][chain$move.type[2:s] == 2]),4), "   rdelt[move 2] accept rate\n")
				if(!missing(groups)) cat(mean(chain$acc.rdelta[2:s][chain$move.type[2:s] == 3]), "   rdelt[move 3] accept rate\n")
				cat(round(colMeans(chain$delta[1:s,ztest ,drop=FALSE]),4), "   post incl probs\n")
				cat(round(colMeans(chain$r[2:s,], na.rm=TRUE),4), "   post mean of r\n")
			}
			print(difftime(Sys.time(), chain$time1))
		}
	}
	chain$time2 <- Sys.time()
	chain$iter <- nsamp
	chain$starting.values <- starting.values
	chain$control.params <- control.params
	chain$X <- X
	chain$Z <- Z
	chain$y <- y
	chain$ztest <- ztest
	chain$data.comps <- data.comps
	if(!missing(Znew)) chain$Znew <- Znew
    class(chain) <- c("bkmrfit", class(chain))
	chain
}

print.bkmrfit <- function(x, q = c(0.025, 0.975), digits = 5, ...) {
    ests <- ExtractEsts(x, q = q)
    summ <- with(ests, rbind(sigsq.eps, beta, r))
    print(round(summ, digits = digits))
}