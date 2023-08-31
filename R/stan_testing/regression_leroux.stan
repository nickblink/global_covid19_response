data {
  int<lower=0> p; // number of variables
  int<lower=0> N;  // number of observations
  int<lower=0> N_F;  // number of facilities
  int<lower=0> N_T;  // number of time points
  matrix[N,p] X; // design matrix
  int y[N];  // output
  vector[p] mu; // prior mean betas
  matrix[p,p] Sigma; // prior variance betas
  matrix[N_F,N_F] W_star; // D - W matrix for Leroux CAR covariance
  matrix[N_F,N_F] I; // Identity matrix
}
parameters {
  real<lower=0> tau2; // CAR variance parameter
  real<lower=0, upper=1> rho; // CAR spatial parameter
  vector[N] phi; // CAR parameter in vector form (not necessary for this implementation but I'll do it anyway).
  vector[p] beta;  // log of rate parameter
}
transformed parameters {
  matrix[N_F,N_F] Q; // Leroux precision matrix
  Q = (1/tau2)*(rho*W_star + (1 - rho)*I); 
}
model {
  beta ~ multi_normal(mu, Sigma); // beta prior
  y ~ poisson_log(X*beta + phi); // likelihood
  for (t in 1:N_T) {
    phi[((t-1)*N_F+1):(t*N_F)] ~ multi_normal_prec(rep_vector(0.0, N_F), Q); // CAR prior
  }
}
generated quantities {
  vector[N] y_exp = exp(X*beta);
  int y_pred[N] = poisson_log_rng(X*beta);
}
