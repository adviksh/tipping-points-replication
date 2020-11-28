functions {
  real normal_density(real x,          // Function argument
                      real xc,         // Complement of function argument
                                       //  on the domain (defined later)
                      real[] theta,    // parameters
                      real[] x_r,      // data (real)
                      int[] x_i) {     // data (integer)
    real xb    = theta[1];
    real sigma = theta[2];
    real delta = theta[3];
    
    real y_obs = x_r[1];
    real x_obs = x_r[2];
    
    return normal_lpdf(y_obs | x_obs < x ? xb : xb + delta, sigma);
  }
}

data {
  int<lower=1> N;
  int<lower=1> K;
  row_vector[K] x_mat[N];
  vector[N] x;
  vector[N] y;
}

transformed data {
  int x_i[0]; // for integrator
}

parameters {
  real alpha;
  vector[K] beta;
  real<lower=40,upper=60> bp;
  real delta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] xb;
  for (n in 1:N) 
    xb[n] = alpha + x_mat[n] * beta;
}

model {
  alpha ~ normal(0, 10);
  beta  ~ normal(0, 10);
  delta ~ normal(0, 10);
  sigma ~ normal(0, 5);
  bp    ~ uniform(40, 60);
  
  for (n in 1:N) {
    target += integrate_1d(normal_density,
                           40,
                           60,
                           { xb[n], sigma, delta},
                           { y[n], x[n] },
                           x_i,
                           0.01);
  }
}

