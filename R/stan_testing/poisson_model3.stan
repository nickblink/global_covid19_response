data {
  int<lower=0> p; // number of variables
  int<lower=0> N;  // number of observations
  matrix[N,p] X; // design matrix
  int<lower=0> y[N];  // data array (counts)
}
parameters {
  vector[p] beta;  // log of rate parameter
}
transformed parameters {
  vector[N] log_lambda;
  log_lambda = X*beta;
}
model {
  y ~ poisson_log(log_lambda);
  // prior
  // beta ~ normal(0, 5);
}
generated quantities {
  vector[N] lambda = exp(log_lambda);
  int yrep[N];
  for (i in 1:N) {
    yrep[i] = poisson_log_rng(log_lambda[i]);
  }
}
