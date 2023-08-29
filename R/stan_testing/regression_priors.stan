data {
  int<lower=0> p; // number of variables
  int<lower=0> N;  // number of observations
  matrix[N,p] X; // design matrix
  int y[N];  // output;
  vector[p] mu; // prior mean betas
  matrix[p,p] Sigma; // prior variance betas
}
parameters {
  vector[p] beta;  // log of rate parameter
}
model {
  beta ~ multi_normal(mu, Sigma); // prior
  y ~ poisson_log(X*beta); // likelihood
}
generated quantities {
  vector[N] y_exp = exp(X*beta);
  int y_pred[N] = poisson_log_rng(X*beta);
}
