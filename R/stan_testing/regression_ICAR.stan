data {
  int<lower=0> p; // number of variables
  int<lower=0> N;  // number of observations
  int<lower=0> F;  // number of facilities
  int<lower=0> T;  // number of time points
  matrix[N,p] X; // design matrix
  int y[N];  // output
  vector[p] mu; // prior mean betas
  matrix[p,p] Sigma; // prior variance betas
  int<lower=0> N_edges; // number of node edges
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
}
parameters {
  vector[N] phi; // CAR parameter in vector form
  vector[p] beta;  // log of rate parameter
}
model {
  beta ~ multi_normal(mu, Sigma); // prior
  y ~ poisson_log(X*beta + phi); // likelihood
  target += -0.5 * dot_self(phi[node1] - phi[node2]); // CAR prior
}
generated quantities {
  vector[N] y_exp = exp(X*beta);
  int y_pred[N] = poisson_log_rng(X*beta);
}
