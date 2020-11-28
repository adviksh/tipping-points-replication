data { 
  int<lower=1> N;  // total number of observations (integer); at least 1
  int<lower=1> K;
  matrix[N,K] x_mat;
  vector[N] x;
  vector[N] y;
}

// the parameters to be estimated from the data
parameters { 
  real alpha;                      // intercept
  vector[K] beta;                  // linear model coefficients
  // real<lower = 0> error;           // standard deviation of residuals
  real<lower = 0, upper = 100> bp; // breakpoint
  real delta;                      // discontinuity at breakpoint
} 

// Functions of estimated parameters.
transformed parameters{
  vector[N] conditional_mean; // the estimated average GJT for each AOA observation
  
  // start with linear prediction
  conditional_mean = alpha + x_mat * beta;
  
  // add delta if above the breakpoint
  for (i in 1:N) 
    conditional_mean[i] = conditional_mean[i] + (x[i] < bp ? 0 : delta);
  
}

// The model itself specifies how the data are expected to have
// been generated and what the prior expectations for the model parameters are.
model {
  alpha ~ normal(0, 50);
  beta ~ normal(0, 100);
  // error ~ normal(0, 5);
  bp ~ normal(15, 0.1);
  delta ~ normal(0, 20);
  y ~ normal(conditional_mean, 5);
}
