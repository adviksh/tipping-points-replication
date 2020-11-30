data { 
  int<lower=1> N; // total number of observations (integer); at least 1
  int<lower=1> K; // polynomial degree
  int<lower=1> B; // number of candidate break points
  
  matrix[N,K] x_mat; // polynomial matrix in x
  
  vector[B] bp; // candidate break points
  vector[N] x; // x values as a vector (for comparing to break points)
  vector[N] y; // y values as a vector
  
  int<lower=1> N_pred;     // number of points to generate predictions for
  vector[N_pred]   x_pred;
  matrix[N_pred,K] x_pred_mat; // x-values to make predictions for
}

transformed data {
  int<lower=0,upper=N> below_bp[B]; // indices of last observation below the break
  for (b in 1:B) {
    below_bp[b] = 1; // initialize to first obsevation
    
    for (n in 2:N) {
      // update if observation n is also below the breakpoint
      below_bp[b] = x[n] < bp[b] ? n : below_bp[b];
    }
  }
}

// the parameters to be estimated from the data
parameters { 
  real alpha;                      // intercept
  vector[K] beta;                  // linear model coefficients
  real delta;                      // discontinuity at breakpoint
  real<lower=0> sigma;             // scale for t distribution
  real<lower=1> nu;                // df for t distribution
} 

// Functions of estimated parameters.
transformed parameters{
  vector[N] xb;
  vector[B] lp;                // log probability for each of the candidates
  xb = alpha + x_mat * beta;
  
  lp = rep_vector(-log(B), B); // implies a uniform prior
  for (b in 1:B) {
    lp[b] = lp[b] + 
      student_t_lpdf(y[1:(below_bp[b])] | nu, xb[1:(below_bp[b])], sigma) +
      student_t_lpdf(y[(below_bp[b] + 1):N] | nu, xb[(below_bp[b] + 1):N] + delta, sigma);
  }
}

model {
  alpha ~ normal(0, 100);
  beta  ~ normal(0, 100);
  delta ~ normal(0, 100);
  sigma ~ normal(0, 100);
  nu ~ gamma(2, 0.1);
  target += log_sum_exp(lp);
}

generated quantities {
  real bp_sim;
  vector[N_pred] y_hat;
  
  bp_sim = bp[categorical_logit_rng(lp)];
  {
    vector[N_pred] mu_pred = alpha + x_pred_mat * beta;
    for (n in 1:N_pred) {
      mu_pred[n] = mu_pred[n] + (x_pred[n] < bp_sim ? 0 : delta);
      y_hat[n]  = student_t_rng(nu, mu_pred[n], sigma);   
    }
  }
}