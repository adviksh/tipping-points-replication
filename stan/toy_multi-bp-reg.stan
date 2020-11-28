data { 
  int<lower=1> N_obs; // total number of observations (integer); at least 1
  int<lower=1> N_cit; // total number of cities
  int<lower=1> N_tract[N_cit]; // number of tracts per city
  
  vector[N_obs] y; // y values as a vector
  
  int<lower=1> K; // polynomial degree
  matrix[N_obs, K] x; // polynomial matrix in x
  
  // Candidate break points
  int<lower=1> B; // number of candidate break points
  vector[B] bp;             
  
   // Indices for the first and last observation below each breakpoint per city
  int<lower=1,upper=N_obs> ct_start_blw[N_cit, B];
  int<lower=1,upper=N_obs> ct_end_blw[N_cit, B];
  
  // Indices for the first and last observation above each breakpoint per city
  int<lower=1,upper=N_obs> ct_start_abv[N_cit, B];
  int<lower=1,upper=N_obs> ct_end_abv[N_cit, B];
  
  // Prior df for student t
  real<lower=0> nu; 
}

transformed data {
  matrix[N_obs, K] Q_ast;
  matrix[K, K] R_ast[N_cit];
  matrix[K, K] R_ast_inverse[N_cit];
  
  for (cc in 1:N_cit) {
    int ct_start;
    int ct_end;
    matrix[N_tract[cc],K] x_ct;
    
    ct_start = ct_start_blw[cc,1];
    ct_end   = ct_end_abv[cc,1];
    x_ct     = x[ct_start:ct_end,];
    
    Q_ast[ct_start:ct_end,] = qr_thin_Q(x_ct) * sqrt(N_tract[cc] - 1);
    R_ast[cc] = qr_thin_R(x_ct) / sqrt(N_tract[cc] - 1);
    R_ast_inverse[cc] = inverse(R_ast[cc]);  
  }
  
}

// the parameters to be estimated from the data
parameters { 
  vector[N_cit] alpha;             // intercept
  vector[K] theta[N_cit];          // linear model coefficients on QR scale
  vector[N_cit] delta;             // discontinuity at breakpoint
  vector<lower=0>[N_cit] sigma;    // scale parameter
} 

// Functions of estimated parameters.
transformed parameters{
  
  vector[B] lp[N_cit];
  {
    for (cc in 1:N_cit) {
      
      vector[N_tract[cc]] xb;
      for (tt in 0:(N_tract[cc] - 1)) {
        xb[tt + 1] = Q_ast[ct_start_blw[cc,1] + tt,] * theta[cc];
      }
      
      // Start with a uniform prior over breakpoints:
      lp[cc,] = rep_vector(-log(B), B);
      
      for (b in 1:B) {
        
        // Possible bug here if all observations fall above, or all fall below
        int N_blw;
        int N_abv;
        N_blw = ct_end_blw[cc,b] - ct_start_blw[cc,b] + 1;
        N_abv = ct_end_abv[cc,b] - ct_start_abv[cc,b] + 1;
        {
          vector[N_blw] xb_blw;
          vector[N_abv] xb_abv;
          xb_blw = xb[1:N_blw];
          xb_abv = xb[(N_blw + 1):(N_blw + N_abv)] + delta[cc];
          lp[cc,b] = lp[cc,b] +
            normal_lpdf(y[ct_start_blw[cc,b]:ct_end_blw[cc,b]] | xb_blw, sigma[cc]) +
            normal_lpdf(y[ct_start_abv[cc,b]:ct_end_abv[cc,b]] | xb_abv, sigma[cc]);
        }
      }
    }  
  }
  
}

model {
  
  alpha ~ normal(0, 50);
  delta ~ normal(0, 20);
  sigma ~ normal(0, 20);
  
  for (cc in 1:N_cit) {
    theta[cc] ~ normal(0, 10);
    target += log_sum_exp(lp[cc]);
  }
  
}

generated quantities {
  int<lower=1,upper=B> bp_idx[N_cit];
  vector[K] beta[N_cit];
  
  for (cc in 1:N_cit) {
    bp_idx[cc] = categorical_logit_rng(lp[cc]);
  }
}