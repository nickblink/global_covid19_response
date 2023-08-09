data {
  int<lower=0> p; // number of variables
  int<lower=0> N;  // number of observations
  matrix[N,p] X; // design matrix
  int y[N];  // output;
}
parameters {
  vector[p] beta;  // log of rate parameter
}
model {
  y ~ poisson_log(X*beta);
}
generated quantities {
  vector[N] y_exp = exp(X*beta);
  int y_pred[N] = poisson_log_rng(X*beta);
}
