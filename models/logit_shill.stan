// models/logit_shill.stan
// Modelo Logit Estándar para clasificación Shill Bidding
// Enlace: p = logistic(X * beta) = inv_logit(eta)

data {
  int<lower=1> N;                     // Número de observaciones
  int<lower=1> P;                     // Número de predictores (incluyendo intercepto)
  matrix[N, P] X;                     // Matriz de diseño estandarizada
  int<lower=0, upper=1> y[N];         // Variable respuesta binaria
  real<lower=0> sigma;                // Desviación estándar a priori para betas
}

parameters {
  vector[P] beta;                     // Coeficientes logísticos
}

model {
  // Prior
  beta ~ normal(0, sigma);

  // Likelihood (Logit estándar) vectorizado con bernoulli_logit para estabilidad numérica
  y ~ bernoulli_logit(X * beta);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N)
    log_lik[i] = bernoulli_logit_lpmf(y[i] | X[i] * beta);
}
