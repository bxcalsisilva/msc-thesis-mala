data {
  int<lower=1> N;              // Observations
  int<lower=1> P;              // Predictors (including intercept)
  matrix[N, P] X;              // Design Matrix
  int<lower=0,upper=1> y[N];   // Binary response
  real<lower=0> sigma;         // Prior SD for betas
}

parameters {
  vector[P] beta;              // Coefficients
  real<lower=0> lambda;        // Power parameter (Renamed from 's' to match paper)
}

model {
  // --- Priors ---
  beta ~ normal(0, sigma);
  
  // Prior for lambda:
  // Matches paper: log(lambda) ~ N(0, 1) -> lambda ~ Lognormal(0, 1) 
  lambda ~ lognormal(0, 1); 

  // --- Likelihood ---
  // Vectorized implementation matches F_p(z) = L(z)^lambda 
  y ~ bernoulli(pow(inv_logit(X * beta), lambda));
}
