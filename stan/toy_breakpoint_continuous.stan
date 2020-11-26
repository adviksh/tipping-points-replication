data {
  int<lower=1> N;
  
  vector[N] y;
  vector[N] x;
}

parameters {
  real mu_lo;
  real mu_hi;
  real sigma;
}

transformed parameters {
  
  real theta[3];
  real z_real[2];
  real z_int[0];
  
  theta[1] = mu_lo;
  theta[2] = mu_hi;
  theta[3] = sigma;
  
  for (ii in 1:N) {
    z_real[1] = y[ii];
    z_real[2] = x[ii];
  }
  
  for (ii in 1:N) {
    vector[N] lp;
    lp[ii] = integrate_1d(normal_density, 0, 100, theta, z_real, z_int);
  }
    
}

model {
  mu_lo ~ normal(0, 2);
  mu_hi ~ normal(0, 2);
  sigma ~ normal(0, 1);
  
  
  target += log_sum_exp(lp);
}

functions {
  real normal_density(real x,          
                      real xc,         
                      real[] theta,    
                      real[] x_r,      
                      int[] x_i) {
  real y_obs = x_r[1];
  real x_obs = x_r[2];
  
  real mu = x_obs > x ? theta[1] : theta[2];
  real sigma = theta[3];

  return 1 / (sqrt(2 * pi()) * sigma) * exp(-0.5 * ((y_obs - mu) / sigma)^2);
}
}