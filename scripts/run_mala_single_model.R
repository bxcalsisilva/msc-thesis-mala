# scripts/run_mala_single_model.R
# Corre UN modelo MALA sobre el dataset Shill Bidding, con las MISMAS metricas
# nuevas que el flujo Stan: LOO-ELPD + Pareto-k (ademas de las de convergencia
# y eficiencia que el driver ya produce).
#
# NOTA sobre Bayes Factor: la verosimilitud marginal es una propiedad del modelo
# (no del sampler). NO se recomputa aqui; el BF se toma como referencia del
# bridge sampling de Stan en run_mala_combine.R. MALA aporta el LOO (propio).
#
# USO:
#   1. Cambia MODEL a "standard", "power" o "reversal"
#   2. Source el script
#   3. Reinicia R (Session > Restart R)
#   4. Repite para los otros dos modelos
#   5. Cuando los tres terminen, corre run_mala_combine.R
#
# Outputs (prefijo 'mala_' / sufijo '_mala_' para NO chocar con los de Stan):
#   output/tables/mala_shill_summary_{model}.csv
#   output/tables/mala_lambda_trace_{model}.csv          (power/reversal)
#   output/loo/loo_mala_{model}.rds
#   output/loo/loo_diagnostics_mala_{model}.csv
#   output/loo/loo_influential_obs_mala_{model}.csv
# ==============================================================================

# ---- CAMBIA ESTO ----
MODEL <- "power"   # "standard" | "power" | "reversal"
# ---------------------

rm(list = setdiff(ls(), "MODEL"))

suppressPackageStartupMessages({
  library(here)
  library(RTMB)
  library(coda)
  library(posterior)
  library(loo)
  library(dplyr)
})

source(here("R", "mala_driver_shill.R"))

stopifnot(MODEL %in% c("standard", "power", "reversal"))

# ==============================================================================
# CONFIGURACION MALA (la establecida en run_shill_thesis_pipeline.R)
# ==============================================================================
PARAMS <- list(
  n_iter   = 440000, # 110000
  n_warmup = 40000, # 10000
  eps_init = 0.0005
)
N_CORES   <- 4          # nº de cadenas (procesos paralelos)
BASE_SEED <- 12345      # reproducibilidad de las semillas por cadena

# Para el LOO: adelgazamos a ~este nº de draws por cadena (memoria controlada).
# log_lik completo (400k draws x 6321 obs) seria ~20 GB; con thinning ~0.4 GB.
THIN_TARGET_POR_CADENA <- 2000

# Covariables FINALES (Fase 1)
COVARIABLES <- c("Bidder_Tendency", "Successive_Outbidding",
                 "Winning_Ratio", "Auction_Duration")

for (d in c("output/tables", "output/loo")) {
  dir.create(here(d), showWarnings = FALSE, recursive = TRUE)
}

# ==============================================================================
# 1. DATOS
# ==============================================================================
cat("Cargando datos...\n")
raw_data <- read.csv(here("data", "Data_shillB.csv"))
X_raw    <- raw_data[, COVARIABLES]
y        <- raw_data$Class
X_scaled <- scale(X_raw)
X_matrix <- cbind(Intercepto = 1, X_scaled)

shill_data <- list(
  X            = X_matrix,
  y            = y,
  p_beta       = ncol(X_matrix),
  sigma2_beta  = 100,    # sd = 10, igual que el prior de Stan
  sigma2_delta = 1,      # no usado con prior uniforme; se conserva por compatibilidad
  param_names  = colnames(X_matrix)
)

cat(sprintf("N=%d | P=%d | positivos=%.1f%%\n\n",
            nrow(X_matrix), ncol(X_matrix), 100 * mean(y == 1)))

# ==============================================================================
# 2. MUESTREO MALA
# ==============================================================================
cat(sprintf("MODELO: %s | cadenas=%d | warmup=%d | iter=%d | eps_init=%.4f\n",
            toupper(MODEL), N_CORES, PARAMS$n_warmup, PARAMS$n_iter, PARAMS$eps_init))

set.seed(BASE_SEED)
res <- run_mala_shill_parallel(
  data_list  = shill_data,
  params     = PARAMS,
  model_type = MODEL,
  n_cores    = N_CORES
)

if (is.null(res)) stop("Muestreo MALA fallido.")

cat(sprintf("Tiempo: %.2f s | Aceptacion: %.1f%% | eps final: %.6f\n",
            res$times$total,
            res$diagnostics$accept_rate * 100,
            res$diagnostics$final_eps))

rhat_max <- max(res$resumen$rhat, na.rm = TRUE)
ess_min  <- min(res$resumen$ess_bulk, na.rm = TRUE)
cat(sprintf("Rhat max: %.4f | ESS min: %.0f\n\n", rhat_max, ess_min))

# ==============================================================================
# 3. GUARDAR ESTIMACIONES DE PARAMETROS
# ==============================================================================
write.csv(
  res$resumen,
  here("output", "tables", sprintf("mala_shill_summary_%s.csv", MODEL)),
  row.names = FALSE
)
cat(sprintf("-> mala_shill_summary_%s.csv\n", MODEL))

# Traza de lambda (power/reversal) — lambda ya viene back-transformada en muestras
if (MODEL != "standard") {
  mcmc <- res$muestras
  lambda_trace <- do.call(rbind, lapply(seq_along(mcmc), function(ch) {
    data.frame(
      iteracion = seq_len(nrow(mcmc[[ch]])),
      cadena    = ch,
      lambda    = as.numeric(mcmc[[ch]][, "lambda"]),
      modelo    = toupper(MODEL)
    )
  }))
  write.csv(
    lambda_trace,
    here("output", "tables", sprintf("mala_lambda_trace_%s.csv", MODEL)),
    row.names = FALSE
  )
  cat(sprintf("-> mala_lambda_trace_%s.csv\n", MODEL))
  rm(lambda_trace)
}

# ==============================================================================
# 4. LOO-ELPD  (reconstruyendo log_lik desde las cadenas de MALA)
# ==============================================================================
cat("\nCalculando log_lik desde las cadenas de MALA (con thinning)...\n")

# 4.1 Extraer draws por cadena, adelgazar y apilar (guardando chain_id)
mcmc <- res$muestras
n_per_chain <- nrow(mcmc[[1]])
thin <- max(1, floor(n_per_chain / THIN_TARGET_POR_CADENA))
idx_thin <- seq(1, n_per_chain, by = thin)

draws_list <- lapply(seq_along(mcmc), function(ch) as.matrix(mcmc[[ch]])[idx_thin, , drop = FALSE])
chain_id   <- rep(seq_along(mcmc), each = length(idx_thin))
all_draws  <- do.call(rbind, draws_list)

p_beta   <- shill_data$p_beta
beta_mat <- all_draws[, 1:p_beta, drop = FALSE]
lambda_v <- if (MODEL == "standard") NULL else all_draws[, "lambda"]

cat(sprintf("  draws usados para LOO: %d (thin=%d, %d por cadena x %d cadenas)\n",
            nrow(all_draws), thin, length(idx_thin), length(mcmc)))

# 4.2 Construir log_lik [S x N] por bloques de draws (memoria controlada)
X <- shill_data$X
yv <- as.integer(shill_data$y)
S  <- nrow(beta_mat)
N  <- nrow(X)
log_lik <- matrix(NA_real_, nrow = S, ncol = N)

chunk <- 1000
clip  <- function(p) pmin(pmax(p, 1e-12), 1 - 1e-12)   # evita log(0)

for (start in seq(1, S, by = chunk)) {
  end  <- min(start + chunk - 1, S)
  rows <- start:end
  eta  <- beta_mat[rows, , drop = FALSE] %*% t(X)        # [chunk x N]

  if (MODEL == "standard") {
    prob <- plogis(eta)
  } else if (MODEL == "power") {
    prob <- plogis(eta) ^ lambda_v[rows]                 # lambda por fila (recicla por columnas)
  } else { # reversal
    prob <- 1 - (plogis(-eta) ^ lambda_v[rows])
  }
  prob <- clip(prob)

  Y <- matrix(yv, nrow = length(rows), ncol = N, byrow = TRUE)
  log_lik[rows, ] <- Y * log(prob) + (1 - Y) * log(1 - prob)
}
rm(eta, prob, Y, all_draws, draws_list, beta_mat); gc()

# 4.3 LOO con r_eff (usa la estructura de cadenas)
r_eff   <- relative_eff(exp(log_lik), chain_id = chain_id, cores = 1)
loo_obj <- loo(log_lik, r_eff = r_eff, cores = 1)
rm(log_lik, r_eff); gc()

print(loo_obj)

k_vals <- pareto_k_values(loo_obj)
cat(sprintf("\nPareto-k (MALA %s):\n", toupper(MODEL)))
cat(sprintf("  k < 0.5      : %d  (%.1f%%)\n", sum(k_vals < 0.5), 100*mean(k_vals < 0.5)))
cat(sprintf("  0.5 - 0.7    : %d  (%.1f%%)\n", sum(k_vals >= 0.5 & k_vals < 0.7), 100*mean(k_vals >= 0.5 & k_vals < 0.7)))
cat(sprintf("  0.7 - 1.0    : %d  (%.1f%%)\n", sum(k_vals >= 0.7 & k_vals < 1.0), 100*mean(k_vals >= 0.7 & k_vals < 1.0)))
cat(sprintf("  k >= 1.0     : %d  (%.1f%%)\n", sum(k_vals >= 1.0), 100*mean(k_vals >= 1.0)))

# ==============================================================================
# 5. GUARDAR LOO
# ==============================================================================
saveRDS(loo_obj, here("output", "loo", sprintf("loo_mala_%s.rds", MODEL)))
cat(sprintf("\n-> loo_mala_%s.rds\n", MODEL))

write.csv(
  data.frame(
    modelo    = toupper(MODEL),
    categoria = c("Bueno (k<0.5)", "Aceptable (0.5-0.7)",
                  "Malo (0.7-1.0)", "Muy malo (k>=1.0)"),
    n_obs     = c(sum(k_vals < 0.5), sum(k_vals >= 0.5 & k_vals < 0.7),
                  sum(k_vals >= 0.7 & k_vals < 1.0), sum(k_vals >= 1.0)),
    pct       = round(c(mean(k_vals < 0.5), mean(k_vals >= 0.5 & k_vals < 0.7),
                        mean(k_vals >= 0.7 & k_vals < 1.0), mean(k_vals >= 1.0)) * 100, 1),
    elpd_loo    = round(loo_obj$estimates["elpd_loo", "Estimate"], 2),
    se_elpd_loo = round(loo_obj$estimates["elpd_loo", "SE"], 2),
    p_loo       = round(loo_obj$estimates["p_loo",    "Estimate"], 2),
    looic       = round(loo_obj$estimates["looic",    "Estimate"], 2)
  ),
  here("output", "loo", sprintf("loo_diagnostics_mala_%s.csv", MODEL)),
  row.names = FALSE
)
cat(sprintf("-> loo_diagnostics_mala_%s.csv\n", MODEL))

idx_inf <- which(k_vals > 0.5)
if (length(idx_inf) > 0) {
  inf_df <- data.frame(obs_idx = idx_inf, k_value = round(k_vals[idx_inf], 4),
                       clase = yv[idx_inf])
  inf_df <- inf_df[order(-inf_df$k_value), ]
  write.csv(inf_df,
            here("output", "loo", sprintf("loo_influential_obs_mala_%s.csv", MODEL)),
            row.names = FALSE)
  cat(sprintf("-> loo_influential_obs_mala_%s.csv  (%d obs con k>0.5)\n", MODEL, nrow(inf_df)))
}

cat(sprintf("\nModelo MALA %s completado.\n", toupper(MODEL)))
cat("Reinicia R (Session > Restart R) y cambia MODEL al siguiente.\n")
cat("Cuando los tres terminen, corre: source(here('scripts','run_mala_combine.R'))\n")
