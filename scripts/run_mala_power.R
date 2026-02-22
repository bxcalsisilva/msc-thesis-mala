# scripts/run_mala_power.R
rm(list = ls())

library(here)
library(progress)
library(dplyr)

source(here("R", "data_generation_power.R")) 
source(here("R", "mala_driver_power.R"))      

# ==============================================================================
# SCENARIO CONFIGURATION (SELECT ONE BLOCK)
# ==============================================================================

# --- A. LOW DIMENSION (10 Params) ---

# [BLOCK 1] Independent
# RHO      <- 0.0; P_TOTAL <- 10

# [BLOCK 2] Medium Correlation
# RHO      <- 0.7; P_TOTAL <- 10

# [BLOCK 3] High Correlation
# RHO      <- 0.9; P_TOTAL <- 10

# --- B. HIGH DIMENSION (50 Params) ---

# [BLOCK 4] Independent
RHO      <- 0.0; P_TOTAL <- 50

# [BLOCK 5] Medium Correlation
# RHO      <- 0.7; P_TOTAL <- 50

# [BLOCK 6] High Correlation
# RHO      <- 0.9; P_TOTAL <- 50

# --- POWER LOGIT SETTINGS ---
LAMBDA_TRUE <- 5

# ==============================================================================

# Configuración del experimento
N_OBS   <- 1000
N_SIM   <- 100
N_CORES <- 4

# Setup de escenario e identificadores
rho_label   <- gsub("^0\\.", "", sprintf("%g", RHO)) 
s_label     <- gsub("^0\\.", "", sprintf("%g", LAMBDA_TRUE)) 
SCENARIO_ID <- sprintf("mala_power_n%d_p%d_rho%s_s%s", N_OBS, P_TOTAL, rho_label, s_label)

cat("Iniciando escenario:", SCENARIO_ID, "\n")

# --- Parámetros reales ---
betas_base <- c(-0.4, 1.5, -1.2, 0.8, -0.9, 0.7, -0.6, 0.25, -0.3, 0.5)
n_noise    <- P_TOTAL - length(betas_base)

set.seed(123) 
betas_noise <- rnorm(n_noise, 0, 0.1)
betas_true  <- c(betas_base, betas_noise)
names(betas_true) <- paste0("beta[", 1:length(betas_true), "]")

# Vector completo para validación (incluye lambda)
betas_true_all <- c(betas_true, "lambda" = LAMBDA_TRUE)

# Configuración MALA
params_mala <- list(
  n_chains = N_CORES, 
  n_warmup = 30000, 
  n_iter   = 330000, 
  eps_init = 0.001 
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
  
  # Generar datos Power Logit
  data_sim <- simulate_power_logit_data(
    n = N_OBS, betas = betas_true, lambda = LAMBDA_TRUE,
    rho = RHO, seed = sim_seeds[i], scale_X = TRUE
  )
  
  mala_data <- list(
    X = data_sim$X, y = data_sim$y, P = length(betas_true),
    sigma2 = 100, param_names = names(betas_true)
  )
  
  # Ejecución
  res <- run_mala_power_logit_parallel(mala_data, params_mala, sim_seeds[i])
  
  if (!is.null(res$resumen)) {
    df <- res$resumen
    
    # Metadatos y diagnósticos
    df$iteration   <- i
    df$scenario    <- SCENARIO_ID
    df$time_total  <- res$times$total
    df$final_eps   <- res$diagnostics$final_eps
    df$accept_rate <- res$diagnostics$accept_rate
    df$multi_ess   <- res$global_metrics$multivariate_ess
    
    # Métricas de precisión
    df$beta_true     <- betas_true_all[match(df$parametro, names(betas_true_all))]
    df$bias          <- df$mean - df$beta_true
    df$squared_error <- df$bias^2
    df$coverage      <- as.integer(df$beta_true >= df$lower_ci & df$beta_true <= df$upper_ci)
    
    results_list[[i]] <- df
  }
  
  if (i %% 10 == 0) gc()
  pb$tick()
}

# --- Guardar resultados ---
final_df <- do.call(rbind, results_list)
write.csv(final_df, sprintf("output/%s_results.csv", SCENARIO_ID), row.names = FALSE)

# Resumen de métricas
report_summary <- final_df %>%
  group_by(parametro) %>%
  summarise(
    Bias = mean(bias),
    RMSE = sqrt(mean(squared_error)),
    Cov  = mean(coverage),
    Rhat = mean(rhat),
    ESS  = mean(ess_bulk),
    Acc  = mean(accept_rate)
  )
print(head(report_summary, 11))