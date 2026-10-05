// models/reversal_power_logit_shill_normal.stan
// Reversal Power Logit para Shill Bidding — VARIANTE CON PRIOR NORMAL en log_lambda.
// Enlace: p = 1 - inv_logit(-eta)^lambda
//
// Diferencia vs reversal_power_logit_shill.stan:
//   - log_lambda SIN limites duros (evita truncamiento del prior)
//   - log_lambda ~ normal(0, 1)  <=>  lambda ~ LogNormal(0, 1)

data {
  int<lower=1> N;                     // Número de observaciones
  int<lower=1> P;                     // Número de predictores (incluyendo intercepto)
  matrix[N, P] X;                     // Matriz de diseño estandarizada
  int<lower=0, upper=1> y[N];         // Variable respuesta binaria
  real<lower=0> sigma;                // Desviación estándar a priori para betas
}

parameters {
  vector[P] beta;                     // Coeficientes logísticos
  real<lower=-10, upper=10> log_lambda; // Limites amplios: solo por seguridad numerica,
                                        // no truncan la normal(0,1) (masa nula mas alla de +-10)
}

transformed parameters {
  real<lower=0> lambda = exp(log_lambda);
}

model {
  // Priores
  beta ~ normal(0, sigma);
  log_lambda ~ normal(0, 1);          // <-- PRIOR NORMAL (antes: uniform(-2, 2))

  // Likelihood (Reversal Power Logit)
  vector[N] eta = X * beta;
  for (i in 1:N) {
    real p = 1 - pow(inv_logit(-eta[i]), lambda);
    y[i] ~ bernoulli(p);
  }
}

generated quantities {
  vector[N] log_lik;
  {
    vector[N] eta = X * beta;
    for (i in 1:N) {
      real p = 1.0 - pow(inv_logit(-eta[i]), lambda);
      log_lik[i] = bernoulli_lpmf(y[i] | p);
    }
  }
}
