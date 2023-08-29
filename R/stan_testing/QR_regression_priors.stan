HOW DO YOU SET A PRIOR ON THE QR DECOMPOSITION?

data {
  int<lower=0> p; // number of variables
  int<lower=0> N;  // number of observations
  matrix[N,p] X; // design matrix
  int y[N];  // output;
  vector u[N]; // prior for betas
}
transformed data {
  matrix[N, p] Q_ast;
  matrix[p, p] R_ast;
  matrix[p, p] R_ast_inverse;
  // thin and scale the QR decomposition
  Q_ast = qr_Q(X)[, 1:p] * sqrt(N - 1);
  R_ast = qr_R(X)[1:p, ] / sqrt(N - 1);
  R_ast_inverse = inverse(R_ast);
}
parameters {
  vector[p] theta;  // coefficients on Q_ast
}
model {
  y ~ poisson_log(Q_ast*theta);
}
generated quantities {
 int y_pred[N];
 vector[p] beta = R_ast_inverse * theta; // is this right? Double check that
 vector[N] y_exp = exp(X*beta);
 y_pred = poisson_log_rng(X*beta);
 //y_pred = poisson_rng(y_exp);
//  int y_pred2[N] = poisson_log_rng(Q_ast*theta);
// vector[N] y_exp2 = exp(Q_ast*theta);
}
