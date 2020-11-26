data {
  int<lower=1> N;
  
  // must be sorted on x
  vector[N] x;
  vector[N] y;
}

transformed data {
  int M;
  real log_unif;
  M = N - 2;
  log_unif = -log(M);
}

parameters {
  real mu_lo;
  real mu_hi;
  real<lower=0>sigma;
}

transformed parameters {
  vector[M] lp;
  {
    vector[M + 1] lp_lo;
    vector[M + 1] lp_hi;
    lp_lo[1] = 0;
    lp_hi[1] = 0;
    for (m in 1:M) {
      lp_lo[m + 1] = lp_lo[m] + normal_lpdf(y[m] | mu_lo, sigma);
      lp_hi[m + 1] = lp_hi[m] + normal_lpdf(y[m] | mu_hi, sigma);
    }
    lp = rep_vector(log_unif + lp_hi[M + 1], M) 
      + head(lp_lo, M) - head(lp_hi, M);
  }
}

model {
  mu_lo ~ normal(0, 2);
  mu_hi ~ normal(0, 2);
  sigma ~ normal(0, 1);
  
  target += log_sum_exp(lp);
}

generated quantities {
  int<lower=1,upper=N> s;
  real delta;
  s = categorical_logit_rng(lp);
  delta = mu_hi - mu_lo;
}