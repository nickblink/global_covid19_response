data {
  int<lower=0> N;  // number of observations
  real y[N];  // data array (counts);
}
parameters {
  real log_lambda;  // log of rate parameter
}
model {
  y ~ poisson_log(log_lambda);
  // prior
  log_lambda ~ normal(0, 5);
}
generated quantities {
  real lambda = exp(log_lambda);
  int yrep[N];
  for (i in 1:N) {
    yrep[i] = poisson_log_rng(log_lambda);
  }
}