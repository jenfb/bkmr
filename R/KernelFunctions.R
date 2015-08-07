GaussianKernel <- function(Z, rho = 1) {
    Kpart <- fields::rdist(Z)^2
    exp(-Kpart/rho)
}

LinearKernel <- function(Z) {
    n <- ncol(Z)
    tcrossprod(Z)/n
}

PolynomialKernel <- function(Z, rho = 1, d = 2) {
    (rho + tcrossprod(Z))^d
}

QuadraticKernel <- function(Z, rho = 1) {
    PolynomialKernel(Z, rho, d = 2)
}
