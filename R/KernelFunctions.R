GaussianKernel <- function(expos, rho = 1) {
    Kpart <- fields::rdist(expos)^2
    exp(-Kpart/rho)
}

LinearKernel <- function(expos) {
    n <- ncol(expos)
    tcrossprod(expos)/n
}

PolynomialKernel <- function(expos, rho = 1, d = 2) {
    (rho + tcrossprod(expos))^d
}

QuadraticKernel <- function(expos, rho = 1) {
    PolynomialKernel(expos, rho, d = 2)
}
