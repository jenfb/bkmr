// Fit a Gaussian process's hyperparameters

data {
  int<lower=1> D; // number of exposure variables
  int<lower=1> K; // number of covariates
  int<lower=1> N; // number of observations
  matrix[N,K] covar;
  vector[D] expos[N];
  vector[N] y; // outcome
}
transformed data {
  vector[N] zero;
  for (i in 1:N)
    zero[i] <- 0;
}
parameters {
  real<lower=0> tau;
  vector<lower=0>[D] r;
  real<lower=0> sigma;
  vector[K] beta; // coefficients for covariates
  vector[N] h;
}
transformed parameters {
  real<lower=0> sigma_sq;
  sigma_sq <- pow(sigma, 2);
}
model {
  matrix[N,N] Sigma;
  vector[N] theta;
  vector[D] rsqrt;

  for(k in 1:D)
    rsqrt[k] <- sqrt(r[k]);

  // off-diagonal elements
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      Sigma[i,j] <- tau * exp(-dot_self((expos[i] - expos[j]) .* rsqrt));
      Sigma[j,i] <- Sigma[i,j];
    }
  }

  // diagonal elements
  for (k in 1:N)
    Sigma[k,k] <- tau;

  r ~ gamma(0.01,0.01);
  tau ~ cauchy(0,1);
  sigma ~ cauchy(0,2.5);

  h ~ multi_normal(zero, Sigma);
  theta <- h + covar * beta;
  y ~ normal(theta, sigma);
}
