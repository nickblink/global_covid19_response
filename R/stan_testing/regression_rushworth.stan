data {
  int<lower=0> p; // number of variables
  int<lower=0> N;  // number of observations
  int<lower=0> N_F;  // number of facilities
  int<lower=0> N_T;  // number of time points
  int<lower=0> N_miss; // number of missing y points
  int<lower=0> N_obs; // number of observed y points
  int<lower=0, upper=N> ind_miss[N_miss]; // indices of missing y points
  int<lower=0, upper=N> ind_obs[N_obs]; // indices of observed y points
  matrix[N,p] X; // design matrix
  int<lower=0> y_obs[N_obs];  // output
  vector[p] mu; // prior mean betas
  matrix[p,p] Sigma; // prior variance betas
  matrix[N_F,N_F] W_star; // D - W matrix for Leroux CAR covariance
  matrix[N_F,N_F] I; // Identity matrix
}
parameters {
  // int<lower=0> y_miss[N_miss];
  real<lower=0> tau2; // CAR variance parameter
  real<lower=0, upper=1> rho; // spatial correlation
  real<lower=0, upper=1> alpha; // temporal correlation
  vector[N] phi; // CAR parameter in vector form (not necessary for this implementation but I'll do it anyway).
  vector[p] beta;  // log of rate parameter
}
transformed parameters {
  matrix[N_F,N_F] Q; // Leroux precision matrix
  Q = (1/tau2)*(rho*W_star + (1 - rho)*I); 
  vector[N_obs] observed_est;
  observed_est = (X*beta)[ind_obs] + phi[ind_obs];
  //int y[N]; // combining the observed and missing data
  //y[ind_obs] = y_obs;
  //y[ind_miss] = y_miss;
}
model {
  beta ~ multi_normal(mu, Sigma); // beta prior
  //y ~ poisson_log(X*beta + phi); // likelihood
  y_obs ~ poisson_log(observed_est);
  phi[1:N_F] ~ multi_normal_prec(rep_vector(0.0, N_F), Q); // MRF initial prior
  for (t in 2:N_T) {
    phi[((t-1)*N_F+1):(t*N_F)] ~ multi_normal_prec(alpha*phi[((t-2)*N_F+1):((t-1)*N_F)], Q); // CAR prior with RF process and Leroux spatial correlation
  }
}
generated quantities {
  vector[N] y_exp = exp(X*beta);
  int y_pred[N] = poisson_log_rng(X*beta);
}
