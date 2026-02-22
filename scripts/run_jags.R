# scripts/run_jags.R
rm(list = ls())

library(here)
library(progress)
library(dplyr)

source(here("R", "data_generation.R"))
source(here("R", "jags_driver.R"))

# ==============================================================================
# SCENARIO CONFIGURATION (SELECT ONE BLOCK)
# ==============================================================================

# --- A. LOW DIMENSION (10 Params) ---

# [BLOCK 1] Independent
# RHO     <- 0.0; P_TOTAL <- 10

# [BLOCK 2] Medium Correlation
# RHO     <- 0.7; P_TOTAL <- 10

# [BLOCK 3] High Correlation
# RHO     <- 0.9; P_TOTAL <- 10

# --- B. HIGH DIMENSION (50 Params) ---

# [BLOCK 4] Independent
# RHO     <- 0.0; P_TOTAL <- 50

# [BLOCK 5] Medium Correlation
# RHO     <- 0.7; P_TOTAL <- 50

# [BLOCK 6] High Correlation
RHO     <- 0.9; P_TOTAL <- 50

# ==============================================================================

# Parámetros globales
N_OBS   <- 500
N_SIM   <- 100
N_CORES <- 4

# Identificador de escenario
rho_label   <- gsub("^0\\.", "", sprintf("%g", RHO))
SCENARIO_ID <- sprintf("jags_n%d_p%d_rho%s", N_OBS, P_TOTAL, rho_label)

cat("Escenario:", SCENARIO_ID, "\n")

# --- Generación de parámetros reales ---
betas_base <- c(-0.4, 1.5, -1.2, 0.8, -0.9, 0.7, -0.6, 0.25, -0.3, 0.5)
n_noise    <- P_TOTAL - length(betas_base)

set.seed(123)
betas_noise <- rnorm(n_noise, 0, 0.1)
betas_true  <- c(betas_base, betas_noise)
names(betas_true) <- paste0("beta[", 1:length(betas_true), "]")

# Configuración JAGS
params_jags <- list(
  n_chains = N_CORES,
  n_adapt  = 500,
  n_burnin = 500,
  n_iter   = 10000
)

# --- Loop de simulación ---
set.seed(123)
sim_seeds <- sample.int(1e6, N_SIM)
results_list <- list()

pb <- progress_bar$new(
  format = " Simulando [:bar] :percent | ETA: :eta",
  total = N_SIM, width = 60
)

for (i in 1:N_SIM) {
  
  # Generar datos y lista para JAGS
  data_sim <- simulate_logit_data(
    n = N_OBS, betas = betas_true, 
    rho = RHO, seed = sim_seeds[i], scale_X = TRUE
  )
  
  jags_data <- list(
    N = N_OBS, P = length(betas_true),
    X = data_sim$X, y = data_sim$y, sigma = 10
  )
  
  # Ejecución
  res <- run_jags_logit_parallel(jags_data, params_jags, sim_seeds[i])
  
  if (!is.null(res$resumen)) {
    df <- res$resumen
    
    # Metadatos y tiempos
    df$iteration  <- i
    df$scenario   <- SCENARIO_ID
    df$time_total <- res$timetotal
    df$multi_ess  <- res$global_metrics$multivariate_ess
    
    # Métricas de precisión y cobertura
    df$beta_true <- betas_true[match(df$parametro, names(betas_true))]
    df$bias      <- df$mean - df$beta_true
    df$sq_error  <- df$bias^2
    df$coverage  <- as.integer(df$beta_true >= df$lower_ci & df$beta_true <= df$upper_ci)
    
    results_list[[i]] <- df
    
    # Guardar cadenas parciales
    saveRDS(res$muestras, 
            file = sprintf("output/chains/jags/%s_iter_%03d.rds", SCENARIO_ID, i), 
            compress = "xz")
  }
  
  if (i %% 10 == 0) gc()
  pb$tick()
}

# --- Exportación final ---
final_df <- do.call(rbind, results_list)
write.csv(final_df, sprintf("output/%s_results.csv", SCENARIO_ID), row.names = FALSE)

cat("\nProceso finalizado. Resultados en:", SCENARIO_ID, "\n")

# Resumen de resultados
report_summary <- final_df %>%
  group_by(parametro) %>%
  summarise(
    Bias = mean(bias),
    RMSE = sqrt(mean(sq_error)),
    Cov  = mean(coverage),
    Rhat = mean(rhat),
    ESS  = mean(ess_bulk)
  )
print(head(report_summary, 10))