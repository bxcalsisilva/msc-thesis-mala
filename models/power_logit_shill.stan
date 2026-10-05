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
  
  // Likelihood (Power Logit)
  vector[N] eta = X * beta;           // Calculamos el predictor lineal de forma vectorizada

  // Likelihood (Power Logit) vectorizado
  // y ~ bernoulli(pow(inv_logit(X * beta), lambda));
  for (i in 1:N) {
    // Aplicamos la asimetría observación por observación
    real p = pow(inv_logit(eta[i]), lambda);
    y[i] ~ bernoulli(p);
  }
}

generated quantities {
  vector[N] log_lik;
  {
    // eta is local to model{}, so we recompute it here inside a local scope
    // The inner { } prevents eta from appearing as a GQ output
    vector[N] eta = X * beta;
    for (i in 1:N) {
      real p = pow(inv_logit(eta[i]), lambda);
      log_lik[i] = bernoulli_lpmf(y[i] | p);
    }
  }
}
