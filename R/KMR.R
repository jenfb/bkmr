## Modified from Maity's code for the function Kernel.reg.Gauss(), available at "http://www4.stat.ncsu.edu/~maity/software/KMTest_garrote_header.r"
## I've changed the bounds for the "lam" parameter (variance parameter), "rho" parameter
## I've added calculations for Var(beta), Var(h) and the function now returns their values
## I also increased the maximum number of iterations
## Finally I changed a few things to increase the speed
#' Kernel machine regression (KMR) with Gaussian kernel
#'
#' The rho parameter is also estimated (rather than fixed)
#'
#' @param y
#' @param expos
#' @param covar
#' @param max.iter
#'
#' @export
FitKMR <- function(y, expos, covar, max.iter = 500, tol = 0.00001) {
## yy: response nx1
## zz: matrix modeled nonparametrically (each col is a sample)
## xx: other variables (each col is a sample)
##
##Estimates: h(z), beta and rho (the kernel parameter)

    yy <- y
    xx <- t(covar)
    zz <- t(expos)

	n = length(yy)

	## iterative algo
	dMat = sqare.diff.mat2(zz)

	bet.hat = matrix(.2, ncol = 1, nrow = nrow(xx))
	########### Step 1(a): fix sigma and rho and estimate lambda
	sig.hat = 1
	rho.hat = 1
	lam.hat = 0.01
	Kmat = kernel.Gauss(dMat, rho.hat)

	conv.crit = tol
	criter = 999
	bet.old = bet.hat
	lam.old = lam.hat
	sig.old = sig.hat
	rho.old = rho.hat
	niter = 0
	#while( (criter > conv.crit) & (niter < 50) ){
	while ((criter > conv.crit) & (niter < max.iter)) {

	#print(c(niter, c(bet.old), lam.old, sig.old, rho.old))
	###########################################
	###### STEP 1 :fix  sigma, lambda and rho and estimate beta and h().
	###########################################
	Kmat = kernel.Gauss(dMat, rho.hat)
	V = sig.hat*diag(1, nrow = n, ncol = n) + lam.hat*Kmat
	#iV = solve(V)
	iV = chol2inv(chol(V))
	#Px = solve(xx%*%iV%*%t(xx))
	Px = chol2inv(chol(xx %*% iV %*% t(xx)))
	#bet.hat = solve(xx%*%iV%*%t(xx)) %*% xx%*%iV%*%yy
	bet.hat = Px %*% xx %*% iV %*% yy
	h.hat = lam.hat*Kmat %*% iV %*% (yy - t(xx) %*% bet.hat)

	###########################################
	###### STEP 2: fix beta and estimate sigma, lambda and rho
	###########################################
	H.mat = (t(xx) %*% Px %*% xx) %*% iV + lam.hat*Kmat %*% iV %*% (diag(1, nrow = n, ncol = n) - t(xx) %*% Px %*% xx %*% iV )

	sig.hat = sum( (yy - t(xx) %*% bet.hat - h.hat)^2 )/ (n - Trace(H.mat))
	if ((n - Trace(H.mat)) < 0)	{
	    sig.hat = (sum((yy - t(xx) %*% bet.hat - h.hat)^2)/ (n))
	}

	lam.hat = optimize(loglik.lam2, interval = c(0.0, 3000), ## I changed this interval to allow larger values ##interval=c(0.0,60),
							Y = (yy - t(xx) %*% bet.hat),
							Kmat = Kmat,
							xx = xx,
							sigmasq = sig.hat)$minimum


	rho.hat = optimize(loglik.rho2, interval = c(0.0, 1000), ## I changed this interval to allow larger values ##interval=c(0.0,40),
							Y = (yy - t(xx) %*% bet.hat),
							dMat = dMat,
							xx = xx,
							lambda = lam.hat,
							sigmasq = sig.hat)$minimum
	##


	criter = sum(abs(bet.hat - bet.old)) +
					abs(lam.hat - lam.old) +
					abs(sig.hat - sig.old) +
					abs(rho.hat - rho.old)

	lik.val = loglik.lam2(lambda = lam.hat, Y = (yy - t(xx) %*% bet.hat), Kmat = Kmat, xx = xx, sigmasq = sig.hat)

	bet.old = bet.hat
	lam.old = lam.hat
	sig.old = sig.hat
	rho.old = rho.hat
	niter = niter + 1


	}
	#print(niter)

	## sanity check
	##if(niter >= 50){
	converge <- TRUE
	if (niter >= max.iter) {
		print("Did not converge")
		converge <- TRUE
		##stop("Did not converge")
	}

	Kmat = kernel.Gauss(dMat, rho.hat)
	V = sig.hat* diag(1, nrow = n, ncol = n) + lam.hat*Kmat
	##iV = solve(V)
	iV <- chol2inv(chol(V))

	##bet.hat = solve(xx%*%iV%*%t(xx)) %*% xx%*%iV%*%yy
	Px = chol2inv(chol(xx %*% iV %*% t(xx)))
	bet.hat = Px %*% xx %*% iV %*% yy
	h.hat = lam.hat*Kmat %*% iV %*% (yy - t(xx) %*% bet.hat)

	## variance-covariance matrix for beta, h.hat, computed using formulas 14,15 of Lin, Liu, and Ghosh
	bet.vcov <- Px
	P <- iV - iV %*% t(xx) %*% Px %*% xx %*% iV
	h.vcov <- lam.hat*Kmat - lam.hat^2*Kmat %*% P %*% Kmat

	## residual df for estimating sigma^2
	H.mat = (t(xx) %*% Px %*% xx)%*%iV +
			lam.hat*Kmat%*%iV%*%( diag(1, nrow=n, ncol=n) - t(xx)%*%Px%*%xx%*%iV )
	df.resid <- ifelse(n - Trace(H.mat) < 0, n, n - Trace(H.mat))

	## variance-covariance matrix of theta=(sigma, tau, rho), computed using formulat given in last line of Section 4.2 of Lin, Liu, and Ghosh
	solve.mat <- iV
	P0.mat <- P
	## this section of code is adapted from the Kernel.test.Gauss function ("http://www4.stat.ncsu.edu/~maity/software/KMTest_garrote_header.r")
	## derivatives of V
	mat.lam = Kmat
	mat.sig = diag(1, nrow = n, ncol = n)
	mat.rho = lam.hat * dMat * Kmat / (rho.hat^2)
	## Info matrix
	I0 = matrix(NA, 3,3)
	I0[1,1] = Trace( P0.mat %*% mat.sig %*% P0.mat %*% mat.sig )/2  ##sig.sig
	I0[1,2] =  Trace( P0.mat %*% mat.sig %*% P0.mat %*% mat.lam )/2  ##sig.lam
	I0[2,1] = I0[1,2]
	I0[1,3] =  Trace( P0.mat %*% mat.sig %*% P0.mat %*% mat.rho )/2  ##sig.rho
	I0[3,1] = I0[1,3]
	I0[2,2] = Trace( P0.mat %*% mat.lam %*% P0.mat %*% mat.lam )/2  ##lam.lam
	I0[2,3] = Trace( P0.mat %*% mat.lam %*% P0.mat %*% mat.rho )/2  ##lam.rho
	I0[3,2] = I0[2,3]
	I0[3,3] =  Trace( P0.mat %*% mat.rho %*% P0.mat %*% mat.rho )/2
	################################
	theta.vcov <- chol2inv(chol(I0))
	rownames(theta.vcov) <- colnames(theta.vcov) <- c("sig","lam","rho")

	return(list(h.hat = h.hat, sig = sig.hat, lam = lam.hat, rho = rho.hat, bet = bet.hat, bet.vcov = bet.vcov, h.vcov = h.vcov, theta.vcov = theta.vcov, df.resid = df.resid, converge = converge))
}

loglik.lam2 = function(lambda, Y, Kmat, xx, sigmasq){

	n = length(Y)
	V = sigmasq* diag(1, nrow = n, ncol = n) + lambda*Kmat
	#iV = solve(V)
	cholV <- chol(V)
	iV = chol2inv(cholV)
	Px <- xx %*% iV %*% t(xx)

	#out = -log(det(V))/2 - log( det(Px) )/2 - t(Y)%*%iV%*%Y/2
	#out = -log(prod(diag(cholV))^2)/2 - log( det(Px) )/2 - t(Y)%*%iV%*%Y/2
	out = -sum(log(diag(cholV))) - log( det(Px) )/2 - t(Y) %*% iV %*% Y/2
	#out = -sum(log(diag(cholV))) - sum(log(diag(chol(Px)))) - t(Y)%*%iV%*%Y/2
	return(-out/n)
}

loglik.rho2 = function(rho, Y, dMat, xx, lambda, sigmasq){

	n = length(Y)
	Kmat = kernel.Gauss(dMat, rho)
	V = sigmasq*diag(1, nrow = n, ncol = n) + lambda*Kmat
	#iV = solve(V)
	cholV <- chol(V)
	iV = chol2inv(cholV)
	Px <- xx %*% iV %*% t(xx)

	#out = -log(det(V))/2 - log( det(xx%*%iV%*%t(xx)) )/2 - t(Y)%*%iV%*%Y/2
	out = -sum(log(diag(cholV))) - log( det(Px) )/2 - t(Y) %*% iV %*% Y/2
	return(-out/n)
}

## function to predict h(Znew) from fit using the above function
predZ <- function(fit, Z, X, y, Znew) {
	if (is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
	if (class(Znew) == "data.frame") Znew <- data.matrix(Znew)
	xx <- t(X)
	yy <- y
	rho <- fit$rho
	sig <- fit$sig
	lam <- fit$lam
	bet <- fit$bet

	n0 <- nrow(Z)
	n1 <- nrow(Znew)
	nall <- n0 + n1
	Kpartall <- as.matrix(dist(rbind(Z, Znew)))^2
	Kmat <- exp(-1/rho*Kpartall)
	Kmat0 <- Kmat[1:n0, 1:n0 ,drop = FALSE]
	Kmat1 <- Kmat[(n0 + 1):nall, (n0 + 1):nall ,drop = FALSE]
	Kmat10 <- Kmat[(n0 + 1):nall, 1:n0 ,drop = FALSE]

	V <- sig*diag(1, nrow = n0, ncol = n0) + lam*Kmat0
	iV <- chol2inv(chol(V))

	hnew <- lam*Kmat10 %*% iV %*% (yy - t(xx) %*% bet)

	## variance-covariance matrix for hnew (using formulas 15 of Lin, Liu, and Ghosh)
	Px <- chol2inv(chol(xx %*% iV %*% t(xx)))
	P <- iV - iV %*% X %*% Px %*% t(X) %*% iV
	hnew.vcov <- lam*Kmat1 - lam^2*Kmat10 %*% P %*% t(Kmat10)

	list(hnew = hnew, hnew.vcov = hnew.vcov)
}

## function to predict h(znew) from fit using the above function
predz <- function(fit, Z, X, y, znew) {
	xx <- t(X)
	yy <- y
	n <- nrow(Z)
	rho <- fit$rho
	sig <- fit$sig
	lam <- fit$lam
	bet <- fit$bet
	Kpart <- as.matrix(dist(Z))^2
	Kmat <- exp(-1/rho*Kpart)

	newKpart <- sapply(1:n, function(x) sum((znew - Z[x,])^2))
	newK <- exp(-1/rho*newKpart)

	V <- sig*diag(1, nrow = n, ncol = n) + lam*Kmat
	iV <- solve(V)

	hnew <- lam*newK %*% iV %*% (yy - t(xx) %*% bet)
	hnew
}

sqare.diff.mat <- function(zz){
    # every column is a sample

    Kmat = matrix(NA, nrow=ncol(zz), ncol=ncol(zz))
    for (ii in 1:ncol(zz)){
        Kmat[ii,] = diag(t(zz-zz[,ii])%*%(zz-zz[,ii]))

    }

    return(Kmat)
}
sqare.diff.mat2 <- function(zz) {
    expos <- t(zz)
    as.matrix(fields::rdist(expos))^2
}

kernel.Gauss <- function(Kmat, scl){

    return(exp(-Kmat/scl))

}
Trace <- function(M){sum(diag(M))}

#' @export
GarroteKernelTest <- function(y, expos, covar, pollutants = 1:ncol(expos), verbose = TRUE, test.method = "pca") {
    df <- data.frame()
    for (i in seq_along(pollutants)) {
        if (verbose) message("testing pollutant ", i)
        df.pol <- GarroteKernelTestPol(y = y, expos = expos, covar = covar, whichpol = i)
        df <- rbind(df, data.frame(pollutant = i, df.pol, converge = attr(df.pol, "converge")))
    }
    df <- dplyr::filter(df, method == test.method) %>%
        dplyr::select(-method)
    attr(df, "method") <- test.method
    df
}
GarroteKernelTestPol <- function(y, expos, covar, whichpol = 1, ...) {
    # test using kernel machine for the effect of the first row of zz

    yy <- y
    xx <- t(covar)
    zz <- t(expos)

    n = length(yy)

    ### estimation under H0
    est.h0 = FitKMR(y = yy, expos = t(zz[-whichpol, ,drop=FALSE]), covar = t(xx), ...)
    converge <- est.h0$converge
    h.hat = est.h0[[1]]
    sig.hat = est.h0[[2]]
    lam.hat = est.h0[[3]]
    rho.hat = est.h0[[4]]
    bet.hat = est.h0[[5]]
    #print(bet.hat)

    dMat = sqare.diff.mat2(zz[-whichpol, ,drop=FALSE])
    Kmat = kernel.Gauss(dMat, rho.hat)

    mat = sig.hat*diag(1, nrow=n, ncol=n) + lam.hat*Kmat
    #solve.mat = solve(mat)
    solve.mat = chol2inv(chol(mat))

    eps.hat = yy - t(xx)%*%bet.hat - h.hat

    ##derivative of the garrote kernel
    mul.mat = matrix(1, ncol=1, nrow=n)%*%zz[whichpol, ]
    mul.mat = -(mul.mat - t(mul.mat))^2
    Kmat.del = mul.mat * Kmat / rho.hat

    ## score stat
    score.chi = lam.hat * t(yy - t(xx)%*%bet.hat) %*% solve.mat %*% Kmat.del %*% solve.mat %*% (yy - t(xx)%*%bet.hat)/2



    ## derivatives of V ('mat' in this code)
    mat.lam = Kmat
    mat.sig = diag(1, nrow=n, ncol=n)
    mat.del = lam.hat*Kmat.del
    mat.rho = lam.hat * dMat * Kmat / (rho.hat^2)

    #Info matrix
    #P0.mat = solve.mat - solve.mat%*%t(xx)%*% solve(xx%*%solve.mat%*%t(xx)) %*% xx%*%solve.mat
    P0.mat = solve.mat - solve.mat%*%t(xx)%*% solve(xx%*%solve.mat%*%t(xx), xx%*%solve.mat)



    #I0 = matrix(NA, 4+nrow(xx),4+nrow(xx))
    I0 = matrix(NA, 4,4)

    I0[1,1] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.del )/2  ##del.del

    I0[1,2] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.sig )/2  ##del.sig
    I0[2,1] = I0[1,2]

    I0[1,3] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.lam )/2  ##del.lam
    I0[3,1] = I0[1,3]

    I0[1,4] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.rho )/2  ##del.rho
    I0[4,1] = I0[1,4]

    #I0[1, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.del%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 1] = I0[1, 5:(4+nrow(xx))]

    I0[2,2] = Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.sig )/2  ##sig.sig

    I0[2,3] =  Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.lam )/2  ##sig.lam
    I0[3,2] = I0[2,3]

    I0[2,4] =  Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.rho )/2  ##sig.rho
    I0[4,2] = I0[2,4]

    #I0[2, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.sig%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 2] = I0[2, 5:(4+nrow(xx))]

    I0[3,3] =  Trace( P0.mat%*%mat.lam%*%P0.mat%*%mat.lam )/2  ##lam.lam

    I0[3,4] =  Trace( P0.mat%*%mat.lam%*%P0.mat%*%mat.rho )/2  ##lam.rho
    I0[4,3] = I0[3,4]

    #I0[3, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.lam%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 3] = I0[3, 5:(4+nrow(xx))]

    I0[4,4] =  Trace( P0.mat%*%mat.rho%*%P0.mat%*%mat.rho )/2

    #I0[4, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.rho%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 4] = I0[4, 5:(4+nrow(xx))]

    #I0[5:(4+nrow(xx)), 5:(4+nrow(xx))] = -xx%*%solve.mat%*%t(xx)


    #Info for delta
    tot.dim = ncol(I0)
    #I.deldel =  I0[1,1] - I0[1,2:tot.dim]%*%solve(I0[2:tot.dim,2:tot.dim])%*%I0[2:tot.dim,1]
    I.deldel =  I0[1,1] - I0[1,2:tot.dim]%*%solve(I0[2:tot.dim,2:tot.dim], I0[2:tot.dim,1])
    #print(I.deldel)

    ####
    ## Chi-square approximation
    md = lam.hat*Trace(Kmat.del%*%P0.mat)/2
    #print(md)

    m.chi = I.deldel / (2*md)
    d.chi = md / m.chi

    pval.chi = 1- pchisq(score.chi/m.chi,d.chi)


    ###################################################
    ############################################
    ### PCA based reconstruction
    ############################################
    ###################################################
    svd.out = eigen(Kmat.del)
    Kmat.del.plus = svd.out$vectors %*% diag(abs(svd.out$values), ncol=length(svd.out$values), nrow=length(svd.out$values)) %*% t(svd.out$vectors)
    #print(svd.out$d)
    #print(eigen(Kmat.del)$values)

    ## derivatives of V ('mat' in this code)
    mat.lam = Kmat
    mat.sig = diag(1, nrow=n, ncol=n)
    mat.del = lam.hat*Kmat.del.plus
    mat.rho = lam.hat * dMat * Kmat / (rho.hat^2)
    ## score stat
    score.pca = t(yy - t(xx)%*%bet.hat) %*% solve.mat %*% mat.del %*% solve.mat %*% (yy - t(xx)%*%bet.hat)/2

    #print(c(score.chi, score.pca))



    #Info matrix
    #P0.mat = solve.mat - solve.mat%*%t(xx)%*% solve(xx%*%solve.mat%*%t(xx)) %*% xx%*%solve.mat
    P0.mat = solve.mat - solve.mat%*%t(xx)%*% solve(xx%*%solve.mat%*%t(xx), xx%*%solve.mat)



    #I0 = matrix(NA, 4+nrow(xx),4+nrow(xx))
    I0 = matrix(NA, 4,4)

    I0[1,1] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.del )/2  ##del.del

    I0[1,2] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.sig )/2  ##del.sig
    I0[2,1] = I0[1,2]

    I0[1,3] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.lam )/2  ##del.lam
    I0[3,1] = I0[1,3]

    I0[1,4] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.rho )/2  ##del.rho
    I0[4,1] = I0[1,4]

    #I0[1, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.del%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 1] = I0[1, 5:(4+nrow(xx))]

    I0[2,2] = Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.sig )/2  ##sig.sig

    I0[2,3] =  Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.lam )/2  ##sig.lam
    I0[3,2] = I0[2,3]

    I0[2,4] =  Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.rho )/2  ##sig.rho
    I0[4,2] = I0[2,4]

    #I0[2, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.sig%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 2] = I0[2, 5:(4+nrow(xx))]

    I0[3,3] =  Trace( P0.mat%*%mat.lam%*%P0.mat%*%mat.lam )/2  ##lam.lam

    I0[3,4] =  Trace( P0.mat%*%mat.lam%*%P0.mat%*%mat.rho )/2  ##lam.rho
    I0[4,3] = I0[3,4]

    #I0[3, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.lam%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 3] = I0[3, 5:(4+nrow(xx))]

    I0[4,4] =  Trace( P0.mat%*%mat.rho%*%P0.mat%*%mat.rho )/2

    #I0[4, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.rho%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 4] = I0[4, 5:(4+nrow(xx))]

    #I0[5:(4+nrow(xx)), 5:(4+nrow(xx))] = -xx%*%solve.mat%*%t(xx)


    #Info for delta
    tot.dim = ncol(I0)
    #I.deldel =  I0[1,1] - I0[1,2:tot.dim]%*%solve(I0[2:tot.dim,2:tot.dim])%*%I0[2:tot.dim,1]
    I.deldel =  I0[1,1] - I0[1,2:tot.dim]%*%solve(I0[2:tot.dim,2:tot.dim], I0[2:tot.dim,1])
    #print(I.deldel)

    ####
    ## Chi-square approximation
    md = Trace(mat.del%*%P0.mat)/2
    #print(md)

    m.pca = I.deldel / (2*md)
    d.pca = md / m.pca

    pval.pca = 1- pchisq(score.pca/m.pca,d.pca)



    ###################################################
    ############################################
    ### PCA based reconstruction (METHOD 2)
    ############################################
    ###################################################
    #svd.out = eigen(Kmat.del)
    Kmat.del.zero = svd.out$vectors %*% diag(svd.out$values*(svd.out$values>0), ncol=length(svd.out$values), nrow=length(svd.out$values)) %*% t(svd.out$vectors)
    #print(svd.out$d)
    #print(eigen(Kmat.del)$values)

    ## derivatives of V ('mat' in this code)
    mat.lam = Kmat
    mat.sig = diag(1, nrow=n, ncol=n)
    mat.del = lam.hat*Kmat.del.zero
    mat.rho = lam.hat * dMat * Kmat / (rho.hat^2)
    ## score stat
    score.pca2 = t(yy - t(xx)%*%bet.hat) %*% solve.mat %*% mat.del %*% solve.mat %*% (yy - t(xx)%*%bet.hat)/2

    #print(c(score.chi, score.pca))



    #Info matrix
    #P0.mat = solve.mat - solve.mat%*%t(xx)%*% solve(xx%*%solve.mat%*%t(xx)) %*% xx%*%solve.mat
    P0.mat = solve.mat - solve.mat%*%t(xx)%*% solve(xx%*%solve.mat%*%t(xx), xx%*%solve.mat)



    #I0 = matrix(NA, 4+nrow(xx),4+nrow(xx))
    I0 = matrix(NA, 4,4)

    I0[1,1] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.del )/2  ##del.del

    I0[1,2] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.sig )/2  ##del.sig
    I0[2,1] = I0[1,2]

    I0[1,3] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.lam )/2  ##del.lam
    I0[3,1] = I0[1,3]

    I0[1,4] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.rho )/2  ##del.rho
    I0[4,1] = I0[1,4]

    #I0[1, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.del%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 1] = I0[1, 5:(4+nrow(xx))]

    I0[2,2] = Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.sig )/2  ##sig.sig

    I0[2,3] =  Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.lam )/2  ##sig.lam
    I0[3,2] = I0[2,3]

    I0[2,4] =  Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.rho )/2  ##sig.rho
    I0[4,2] = I0[2,4]

    #I0[2, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.sig%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 2] = I0[2, 5:(4+nrow(xx))]

    I0[3,3] =  Trace( P0.mat%*%mat.lam%*%P0.mat%*%mat.lam )/2  ##lam.lam

    I0[3,4] =  Trace( P0.mat%*%mat.lam%*%P0.mat%*%mat.rho )/2  ##lam.rho
    I0[4,3] = I0[3,4]

    #I0[3, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.lam%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 3] = I0[3, 5:(4+nrow(xx))]

    I0[4,4] =  Trace( P0.mat%*%mat.rho%*%P0.mat%*%mat.rho )/2

    #I0[4, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.rho%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 4] = I0[4, 5:(4+nrow(xx))]

    #I0[5:(4+nrow(xx)), 5:(4+nrow(xx))] = -xx%*%solve.mat%*%t(xx)


    #Info for delta
    tot.dim = ncol(I0)
    #I.deldel =  I0[1,1] - I0[1,2:tot.dim]%*%solve(I0[2:tot.dim,2:tot.dim])%*%I0[2:tot.dim,1]
    I.deldel =  I0[1,1] - I0[1,2:tot.dim]%*%solve(I0[2:tot.dim,2:tot.dim], I0[2:tot.dim,1])
    #print(I.deldel)

    ####
    ## Chi-square approximation
    md = Trace(mat.del%*%P0.mat)/2
    #print(md)

    m.pca2 = I.deldel / (2*md)
    d.pca2 = md / m.pca2

    pval.pca2 = 1- pchisq(score.pca2/m.pca2, d.pca2)

    ret <- data.frame(method = c("chisq", "pca", "pca2"),
                      test.stat = c(score.chi, score.pca, score.pca2),
                      p.val = c(pval.chi, pval.pca, pval.pca2),
                      m = c(m.chi, m.pca, m.pca2),
                      d = c(d.chi, d.pca, d.pca2),
                      stringsAsFactors = FALSE
                      )
    attr(ret, "converge") <- converge
    ret
}
