data {
  int<lower=1> N;              // Número de observaciones
  int<lower=1> P;              // Número de predictores (incluye intercepto)
  matrix[N, P] X;              // Matriz de diseño
  int<lower=0,upper=1> y[N];   // Variable de respuesta binaria
  real<lower=0> sigma;         // Desviación estándar del prior
}
parameters {
  vector[P] beta;              // Vector de coeficientes
}
model {
  // Prior: beta ~ Normal(0, sigma)
  beta ~ normal(0, sigma);

  // Verosimilitud (Likelihood)
  y ~ bernoulli_logit(X * beta);
}
