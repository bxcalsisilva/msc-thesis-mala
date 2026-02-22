# 1. Load functions
library(here)
source(here("R", "data_generation.R"))
source(here("R", "jags_driver.R"))

# 2. Generate test data
set.seed(123)
betas_test <- c(-0.5, 0.8, -0.2) # Intercept + 2 predictors
names(betas_test) <- paste0("beta[", 1:3, "]")

test_data <- simulate_logit_data(n = 1000, betas = betas_test, rho = 0, seed = 123)

# 3. Frequentist Fit (Standard R GLM)
# We use -1 because 'test_data$X' already contains the column of 1s for the intercept
fit_glm <- glm(test_data$y ~ test_data$X - 1, family = binomial()) 

# 4. Bayesian Fit (JAGS)
jags_data <- list(N=1000, P=3, X=test_data$X, y=test_data$y, sigma=10)
params <- list(n_chains=4, n_adapt=500, n_burnin=500, n_iter=10000)

fit_jags <- run_jags_logit_parallel(jags_data, params, base_seed=123)

# 5. Final Comparison
if(!is.null(fit_jags$resumen)) {
  comparison <- data.frame(
    Parameter  = names(betas_test),
    True_Value = betas_test,
    GLM_Est    = coef(fit_glm),
    JAGS_Mean  = fit_jags$resumen$mean,
    Difference = abs(coef(fit_glm) - fit_jags$resumen$mean)
  )
  
  print(comparison)
} else {
  cat("Error: JAGS returned no results")
}