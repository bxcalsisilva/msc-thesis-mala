# scripts/run_stan_single_model.R
# Corre UN modelo Stan sobre el dataset Shill Bidding.
#
# USO:
#   1. Cambia MODEL a "standard", "power" o "reversal"
#   2. Source el script (o Ctrl+Shift+Enter en RStudio)
#   3. Reinicia R (Session > Restart R)
#   4. Repite para los otros dos modelos
#   5. Cuando los tres terminen, corre run_stan_combine.R
#
# Cada ejecucion escribe en output/loo/ y output/tables/ sin tocar
# los resultados de los otros modelos.
# ==============================================================================

# ---- CAMBIA ESTO ----
MODEL <- "reversal"   # "standard" | "power" | "reversal"
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

stopifnot(MODEL %in% c("standard", "power", "reversal"))

# ==============================================================================
# PARAMETROS MCMC — ajustados por modelo para controlar uso de RAM
#
# Memoria pico (extract_log_lik crea 2 copias de log_lik):
#   standard: 4 cadenas × 4000 draws × 6321 obs × 8 bytes × 2 = ~1.6 GB
#   power   : 2 cadenas × 4000 draws × 6321 obs × 8 bytes × 2 = ~0.8 GB
#   reversal: 2 cadenas × 4000 draws × 6321 obs × 8 bytes × 2 = ~0.8 GB
# ==============================================================================
PARAMS_BY_MODEL <- list(
  standard = list(n_chains = 4, n_warmup = 1000, n_iter = 5000,
                  adapt_delta = 0.90, max_treedepth = 12),
  power    = list(n_chains = 2, n_warmup = 1000, n_iter = 6000,
                  adapt_delta = 0.99, max_treedepth = 14),
  reversal = list(n_chains = 2, n_warmup = 1000, n_iter = 7000,
                  adapt_delta = 0.99, max_treedepth = 14)
)

PARAMS    <- PARAMS_BY_MODEL[[MODEL]]
BASE_SEED <- 12345

stan_files <- list(
  standard = here("models", "logit_shill.stan"),
  power    = here("models", "power_logit_shill.stan"),
  reversal = here("models", "reversal_power_logit_shill.stan")
)

for (d in c("output/tables", "output/loo")) {
  dir.create(here(d), showWarnings = FALSE, recursive = TRUE)
}

# ==============================================================================
# 1. DATOS
# ==============================================================================
cat("Cargando datos...\n")
raw_data <- read.csv(here("data", "Data_shillB.csv"))
# Covariables FINALES (Fase 1: seleccion por IC 95% + corroboracion LOO).
# Para volver al conjunto completo, usa: X_raw <- raw_data[, 4:12]
COVARIABLES <- c("Bidder_Tendency", "Successive_Outbidding",
                 "Winning_Ratio", "Auction_Duration")
X_raw    <- raw_data[, COVARIABLES]
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
  here("output", "tables", sprintf("stan_shill_summary_%s.csv", MODEL)),
  row.names = FALSE
)
cat(sprintf("-> stan_shill_summary_%s.csv\n", MODEL))

# Para Power y Reversal: guardar traza de lambda como CSV
if (MODEL != "standard") {
  draws_df <- as_draws_df(as_draws(res$muestras))
  lambda_trace <- data.frame(
    iteracion = draws_df$.iteration,
    cadena    = as.integer(draws_df$.chain),
    lambda    = as.numeric(draws_df$lambda),
    modelo    = toupper(MODEL)
  )
  write.csv(
    lambda_trace,
    here("output", "tables", sprintf("stan_lambda_trace_%s.csv", MODEL)),
    row.names = FALSE
  )
  cat(sprintf("-> stan_lambda_trace_%s.csv\n", MODEL))
  rm(draws_df, lambda_trace)
  gc()
}

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
  saveRDS(bs_obj, here("output", "loo", sprintf("bridge_%s.rds", MODEL)))
  cat(sprintf("-> bridge_%s.rds\n", MODEL))
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
saveRDS(loo_obj, here("output", "loo", sprintf("loo_%s.rds", MODEL)))
cat(sprintf("\n-> loo_%s.rds\n", MODEL))

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
  here("output", "loo", sprintf("loo_diagnostics_%s.csv", MODEL)),
  row.names = FALSE
)
cat(sprintf("-> loo_diagnostics_%s.csv\n", MODEL))

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
    here("output", "loo", sprintf("loo_influential_obs_%s.csv", MODEL)),
    row.names = FALSE
  )
  cat(sprintf("-> loo_influential_obs_%s.csv  (%d obs con k>0.5)\n",
              MODEL, nrow(inf_df)))
}

cat(sprintf("\nModelo %s completado.\n", toupper(MODEL)))
cat("Ahora reinicia R (Session > Restart R) y cambia MODEL al siguiente.\n")
cat("Cuando los tres terminen, corre: source(here('scripts','run_stan_combine.R'))\n")
