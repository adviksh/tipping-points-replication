data {
  int<lower=0> N;   // number of observations
  int<lower=0> K;   // number of predictors
  real<lower=0> nu; // prior df
  matrix[N, K] x;   // predictors
  vector[N] y;      // outcome
}

transformed data {
  matrix[N, K] Q_ast;
  matrix[K, K] R_ast;
  matrix[K, K] R_ast_inverse;
  Q_ast = qr_thin_Q(x) * sqrt(N - 1);
  R_ast = qr_thin_R(x) / sqrt(N - 1);
  R_ast_inverse = inverse(R_ast);
}

parameters {
  real alpha;          // intercept
  vector[K] theta;     // coefficients on Q_ast
  real<lower=0> sigma; // error scale
}

model {
  y ~ student_t(nu, alpha + Q_ast * theta, sigma);
  theta ~ normal(0, 1);
  sigma ~ normal(0, 2);
}

generated quantities {
  vector[K] beta;
  beta = R_ast_inverse * theta;
}
