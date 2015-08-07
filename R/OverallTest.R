#' Overall kernel test
#'
#' Conduct an overall test of the mixture using the approach of Liu, Lin, and Ghosh (2007) Biometrics.
#'
#' @param rho tuning parameter for various kernel functions
#' @param y
#' @param Z
#' @param X
#' @param kernel
#' @param d tuning parameter for dth polynomial kernel
#' @export
OverallTest <- function(y, Z, X = NULL, kernel = "gaussian", rho = 1, d = 2) {
    
    n = length(y)
    p = ncol(X)
    m = ncol(Z)

    if (!is.null(X)) {
        X <- cbind(1, X)
    } else {
        X <- as.matrix(rep(1, n))
    }

    if (inherits(kernel, "matrix")) {
        K <- kernel
    } else {
        if (tolower(kernel) == "gaussian") {
            K <- GaussianKernel(Z = Z, rho = rho)
        } else if (kernel == "linear") {
            K <- LinearKernel(Z = Z)
        } else if (kernel == "quadratic") {
            K <- QuadraticKernel(Z = Z)
        } else if (kernel == "polynomal") {
            K <- PolynomialKernel(Z = Z, rho = rho, d = d)
        }
    }

    mod <- lm(y ~ 0 + X)
    sig2 <- summary(mod)$sigma^2
    res <- resid(mod)

    Qt <- drop(t(res) %*% K %*% res)/2/sig2
    P0 <- diag(n) - X %*% solve(crossprod(X), t(X))
    P0K <- P0 %*% K

    Itt <- sum(diag(P0K %*% P0K))/2
    Its <- sum(diag(P0K %*% P0))/2
    Iss <- sum(diag(P0 %*% P0))/2
    etilde <- sum(diag(P0K))/2
    Ihat <- Itt - Its^2/Iss
    kappa <- Ihat/2/etilde
    nu <- 2*etilde^2/Ihat

    S <- Qt/kappa
    p <- pchisq(S, df = nu, lower.tail = FALSE)

    ret <- list(chisq = drop(S), df = nu, Q = drop(Qt), scale = kappa, p.val = p)
    ret
}

