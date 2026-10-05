# scripts/run_stan_single_model_normal.R
# EXPERIMENTO PRIOR NORMAL — corre UN modelo asimetrico Stan sobre Shill Bidding
# usando log_lambda ~ normal(0,1) (sin limites duros), en lugar de uniform(-2,2).
#
# Es una copia independiente de run_stan_single_model.R:
#   - NO toca el flujo vigente ni sus outputs.
#   - Solo aplica a "power" y "reversal" (el "standard" no tiene lambda; su
#     resultado es invariante al prior de lambda y se reutiliza desde el flujo
#     vigente al combinar con run_stan_combine_normal.R).
#   - Todos los outputs llevan sufijo "_normal":
#       output/tables/stan_shill_summary_{power,reversal}_normal.csv
#       output/tables/stan_lambda_trace_{power,reversal}_normal.csv
#       output/loo/{loo,bridge}_{power,reversal}_normal.rds
#       output/loo/loo_diagnostics_{power,reversal}_normal.csv
#       output/loo/loo_influential_obs_{power,reversal}_normal.csv
#
# USO:
#   1. Cambia MODEL a "power" o "reversal"
#   2. Source el script
#   3. Reinicia R (Session > Restart R)
#   4. Repite para el otro modelo
#   5. Cuando ambos terminen, corre run_stan_combine_normal.R
# ==============================================================================

# ---- CAMBIA ESTO ----
MODEL <- "reversal"   # "power" | "reversal"
# ---------------------

rm(list = setdiff(ls(), "MODEL"))

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(posterior)
  library(loo)
  library(dplyr)
  library(bridgesampling)
})

source(here("R", "stan_driver_shill.R"))
rstan_options(auto_write = TRUE)

stopifnot(MODEL %in% c("power", "reversal"))

# Sufijo de salida: separa por completo de los resultados con prior uniforme.
TAG <- paste0(MODEL, "_normal")

# ==============================================================================
# PARAMETROS MCMC (identicos al flujo vigente para que la comparacion sea justa)
# ==============================================================================
PARAMS_BY_MODEL <- list(
  power    = list(n_chains = 2, n_warmup = 1000, n_iter = 5000,
                  adapt_delta = 0.99, max_treedepth = 14),
  reversal = list(n_chains = 2, n_warmup = 1000, n_iter = 5000,
                  adapt_delta = 0.99, max_treedepth = 14)
)

PARAMS    <- PARAMS_BY_MODEL[[MODEL]]
BASE_SEED <- 12345

# Modelos Stan con PRIOR NORMAL (sin limites duros en log_lambda)
stan_files <- list(
  power    = here("models", "power_logit_shill_normal.stan"),
  reversal = here("models", "reversal_power_logit_shill_normal.stan")
)

cat(sprintf("EXPERIMENTO PRIOR NORMAL | modelo=%s | prior log_lambda ~ normal(0,1)\n",
            toupper(MODEL)))
cat(sprintf("Outputs -> *_%s.*\n", TAG))

for (d in c("output/tables", "output/loo")) {
  dir.create(here(d), showWarnings = FALSE, recursive = TRUE)
}

# ==============================================================================
# 1. DATOS
# ==============================================================================
cat("Cargando datos...\n")
raw_data <- read.csv(here("data", "Data_shillB.csv"))
X_raw    <- raw_data[, 4:12]
y        <- raw_data$Class
X_scaled <- scale(X_raw)
X_matrix <- cbind(Intercepto = 1, X_scaled)

stan_data <- list(
  N     = nrow(X_matrix),
  P     = ncol(X_matrix),
  X     = X_matrix,
  y     = as.integer(y),
  sigma = 10,
  param_names = colnames(X_matrix)
)

cat(sprintf("N=%d | P=%d | positivos=%.1f%%\n\n",
            stan_data$N, stan_data$P, 100 * mean(stan_data$y == 1)))

# ==============================================================================
# 2. COMPILAR
# ==============================================================================
cat(sprintf("Compilando %s...\n", basename(stan_files[[MODEL]])))
compiled_model <- rstan::stan_model(file = stan_files[[MODEL]])

# ==============================================================================
# 3. MUESTREO
# ==============================================================================
cat(sprintf("\nMODELO: %s | cadenas=%d | warmup=%d | iter=%d\n",
            toupper(MODEL), PARAMS$n_chains, PARAMS$n_warmup, PARAMS$n_iter))

res <- run_stan_shill_parallel(
  data_list  = stan_data,
  params     = PARAMS,
  model_obj  = compiled_model,
  base_seed  = BASE_SEED,
  model_type = MODEL
)

if (is.null(res)) stop("Muestreo fallido.")

# Liberar modelo compilado — ya no se necesita
rm(compiled_model)
gc()

cat(sprintf("Tiempo: %.2f s | Aceptacion: %.1f%% | StepSize: %.5f\n",
            res$times$total,
            res$diagnostics$accept_rate * 100,
            res$diagnostics$final_stepsize))

rhat_max <- max(res$resumen$rhat, na.rm = TRUE)
ess_min  <- min(res$resumen$ess_bulk, na.rm = TRUE)
cat(sprintf("Rhat max: %.4f | ESS min: %.0f\n\n", rhat_max, ess_min))

# ==============================================================================
# 4. GUARDAR ESTIMACIONES DE PARAMETROS
# ==============================================================================
write.csv(
  res$resumen,
  here("output", "tables", sprintf("stan_shill_summary_%s.csv", TAG)),
  row.names = FALSE
)
cat(sprintf("-> stan_shill_summary_%s.csv\n", TAG))

# Traza de lambda como CSV (clave: aqui se ve si lambda se despega de la frontera)
draws_df <- as_draws_df(as_draws(res$muestras))
lambda_trace <- data.frame(
  iteracion = draws_df$.iteration,
  cadena    = as.integer(draws_df$.chain),
  lambda    = as.numeric(draws_df$lambda),
  modelo    = toupper(MODEL)
)
write.csv(
  lambda_trace,
  here("output", "tables", sprintf("stan_lambda_trace_%s.csv", TAG)),
  row.names = FALSE
)
cat(sprintf("-> stan_lambda_trace_%s.csv\n", TAG))
rm(draws_df, lambda_trace)
gc()

# ==============================================================================
# 4b. BRIDGE SAMPLING — log-verosimilitud marginal para Bayes Factor
# ==============================================================================
cat("\nCalculando Bridge Sampling (log-verosimilitud marginal)...\n")

bs_obj <- tryCatch({
  bridge_sampler(res$muestras, silent = FALSE)
}, error = function(e) {
  message("  Bridge sampling fallido: ", e$message)
  return(NULL)
})

if (!is.null(bs_obj)) {
  logml_val <- bs_obj$logml
  err_val   <- tryCatch(error_measures(bs_obj)$re2, error = function(e) NA)
  cat(sprintf("  log(Z) = %.4f", logml_val))
  if (!is.na(err_val)) cat(sprintf("  |  RE^2 = %.6f", err_val))
  cat("\n")
  saveRDS(bs_obj, here("output", "loo", sprintf("bridge_%s.rds", TAG)))
  cat(sprintf("-> bridge_%s.rds\n", TAG))
} else {
  cat("  Bridge sampling no disponible para este modelo.\n")
}

# ==============================================================================
# 5. LOO-ELPD
# ==============================================================================
cat("\nCalculando LOO-ELPD...\n")

log_lik <- extract_log_lik(
  res$muestras,
  parameter_name = "log_lik",
  merge_chains   = FALSE
)

# Liberar el fit — ya no se necesita, solo necesitamos log_lik
res$muestras <- NULL
gc()

rel_eff <- relative_eff(exp(log_lik), cores = 1)
loo_obj <- loo(log_lik, r_eff = rel_eff, cores = 1)

rm(log_lik, rel_eff)
gc()

print(loo_obj)

k_vals <- pareto_k_values(loo_obj)
cat(sprintf("\nPareto-k (%s):\n", toupper(MODEL)))
cat(sprintf("  k < 0.5      : %d  (%.1f%%)\n",
            sum(k_vals < 0.5),                       100 * mean(k_vals < 0.5)))
cat(sprintf("  0.5 - 0.7    : %d  (%.1f%%)\n",
            sum(k_vals >= 0.5 & k_vals < 0.7),       100 * mean(k_vals >= 0.5 & k_vals < 0.7)))
cat(sprintf("  0.7 - 1.0    : %d  (%.1f%%)\n",
            sum(k_vals >= 0.7 & k_vals < 1.0),       100 * mean(k_vals >= 0.7 & k_vals < 1.0)))
cat(sprintf("  k >= 1.0     : %d  (%.1f%%)\n",
            sum(k_vals >= 1.0),                      100 * mean(k_vals >= 1.0)))

# ==============================================================================
# 6. GUARDAR LOO (RDS para combine, CSV para lectura directa)
# ==============================================================================
saveRDS(loo_obj, here("output", "loo", sprintf("loo_%s.rds", TAG)))
cat(sprintf("\n-> loo_%s.rds\n", TAG))

# Diagnosticos Pareto-k resumidos
write.csv(
  data.frame(
    modelo    = toupper(MODEL),
    categoria = c("Bueno (k<0.5)", "Aceptable (0.5-0.7)",
                  "Malo (0.7-1.0)", "Muy malo (k>=1.0)"),
    n_obs     = c(sum(k_vals < 0.5),
                  sum(k_vals >= 0.5 & k_vals < 0.7),
                  sum(k_vals >= 0.7 & k_vals < 1.0),
                  sum(k_vals >= 1.0)),
    pct       = round(c(mean(k_vals < 0.5),
                        mean(k_vals >= 0.5 & k_vals < 0.7),
                        mean(k_vals >= 0.7 & k_vals < 1.0),
                        mean(k_vals >= 1.0)) * 100, 1),
    elpd_loo    = round(loo_obj$estimates["elpd_loo", "Estimate"], 2),
    se_elpd_loo = round(loo_obj$estimates["elpd_loo", "SE"], 2),
    p_loo       = round(loo_obj$estimates["p_loo",    "Estimate"], 2),
    looic       = round(loo_obj$estimates["looic",    "Estimate"], 2)
  ),
  here("output", "loo", sprintf("loo_diagnostics_%s.csv", TAG)),
  row.names = FALSE
)
cat(sprintf("-> loo_diagnostics_%s.csv\n", TAG))

# Observaciones influyentes (k > 0.5)
idx_inf <- which(k_vals > 0.5)
if (length(idx_inf) > 0) {
  inf_df <- data.frame(
    obs_idx = idx_inf,
    k_value = round(k_vals[idx_inf], 4),
    clase   = stan_data$y[idx_inf]
  )
  inf_df <- inf_df[order(-inf_df$k_value), ]
  write.csv(
    inf_df,
    here("output", "loo", sprintf("loo_influential_obs_%s.csv", TAG)),
    row.names = FALSE
  )
  cat(sprintf("-> loo_influential_obs_%s.csv  (%d obs con k>0.5)\n",
              TAG, nrow(inf_df)))
}

cat(sprintf("\nModelo %s (prior normal) completado.\n", toupper(MODEL)))
cat("Ahora reinicia R (Session > Restart R) y cambia MODEL al otro modelo.\n")
cat("Cuando ambos terminen, corre: source(here('scripts','run_stan_combine_normal.R'))\n")
