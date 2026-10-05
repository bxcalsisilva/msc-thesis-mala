// models/reversal_power_logit_shill.stan
// Modelo Reversal Power Logit para clasificación Shill Bidding
// Enlace: p = 1 - logistic(-eta)^lambda = 1 - inv_logit(-eta)^lambda

data {
  int<lower=1> N;                     // Número de observaciones
  int<lower=1> P;                     // Número de predictores (incluyendo intercepto)
  matrix[N, P] X;                     // Matriz de diseño estandarizada
  int<lower=0, upper=1> y[N];         // Variable respuesta binaria
  real<lower=0> sigma;                // Desviación estándar a priori para betas
}

parameters {
  vector[P] beta;                     // Coeficientes logísticos
  real<lower=-2, upper=2> log_lambda; // Prior uniforme en este rango
}

transformed parameters {
  real<lower=0> lambda = exp(log_lambda);
}

model {
  // Priores
  beta ~ normal(0, sigma);

  // La uniforme es implícita por los límites; la línea siguiente es opcional:
  log_lambda ~ uniform(-2, 2);

  // Likelihood (Reversal Power Logit)
  vector[N] eta = X * beta;

  for (i in 1:N) {
    // p = 1 - logistic(-eta)^lambda
    real p = 1 - pow(inv_logit(-eta[i]), lambda);
    y[i] ~ bernoulli(p);
  }
}

generated quantities {
  vector[N] log_lik;
  {
    // eta is local to model{}, so we recompute it here inside a local scope
    vector[N] eta = X * beta;
    for (i in 1:N) {
      real p = 1.0 - pow(inv_logit(-eta[i]), lambda);
      log_lik[i] = bernoulli_lpmf(y[i] | p);
    }
  }
}
