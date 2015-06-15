GaussianKernel <- function(expos, rho = 1) {
    Kpart <- fields::rdist(expos)^2
    exp(-Kpart/rho)
}

LinearKernel <- function(expos) {
    n <- ncol(expos)
    tcrossprod(expos)/n
}

PolynomialKernel <- function(expos, rho, d) {
    (rho + tcrossprod(expos))^d
}

QuadraticKernel <- function(expos, rho) {
    PolynomialKernel(expos, rho, d = 2)
}
