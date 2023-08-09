data {
  int<lower=0> p; // number of variables
  int<lower=0> N;  // number of observations
  matrix[N,p] X; // design matrix
  int y[N];  // output;
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
  vector[p] beta = R_ast_inverse * theta;
  vector[N] y_exp = exp(X*beta);
  int y_pred[N] = poisson_log_rng(X*beta);
  int y_pred2[N] = poisson_log_rng(Q_ast*theta);
  vector[N] y_exp2 = exp(Q_ast*theta);
}
