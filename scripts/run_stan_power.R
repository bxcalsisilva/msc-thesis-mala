# scripts/run_stan_power.R
rm(list = ls())

library(here)
library(progress)
library(dplyr)

source(here("R", "data_generation_power.R"))
source(here("R", "stan_driver_power.R"))

# ==============================================================================
# 1. SCENARIO CONFIGURATION (SELECT ONE BLOCK)
# ==============================================================================

# --- A. LOW DIMENSION (10 Params) ---

# [BLOCK 1] Independent
# RHO     <- 0.0; P_TOTAL <- 10

# [BLOCK 2] Medium Correlation
# RHO     <- 0.7; P_TOTAL <- 10

# [BLOCK 3] High Correlation (ACTIVE)
# RHO     <- 0.9; P_TOTAL <- 10

# --- B. HIGH DIMENSION (50 Params) ---

# [BLOCK 4] Independent
# RHO     <- 0.0; P_TOTAL <- 50

# [BLOCK 5] Medium Correlation
# RHO     <- 0.7; P_TOTAL <- 50

# [BLOCK 6] High Correlation
RHO     <- 0.9; P_TOTAL <- 50

# --- POWER LOGIT SPECIFIC ---
LAMBDA_TRUE  <- 0.8   

# ==============================================================================

# Configuración del experimento
N_OBS   <- 1000
N_SIM   <- 100
N_CORES <- 4     

# Identificadores de escenario
rho_label   <- gsub("^0\\.", "", sprintf("%g", RHO)) 
s_label     <- gsub("^0\\.", "", sprintf("%g", LAMBDA_TRUE)) 
SCENARIO_ID <- sprintf("stan_power_n%d_p%d_rho%s_s%s", N_OBS, P_TOTAL, rho_label, s_label)

cat("Ejecutando escenario:", SCENARIO_ID, "\n")

# --- Generación de parámetros reales ---
betas_base <- c(-0.4, 1.5, -1.2, 0.8, -0.9, 0.7, -0.6, 0.25, -0.3, 0.5)
n_noise    <- P_TOTAL - length(betas_base)

set.seed(123) 
betas_noise <- if(n_noise > 0) rnorm(n_noise, 0, 0.1) else NULL
betas_true  <- c(betas_base, betas_noise)
names(betas_true) <- paste0("beta[", 1:length(betas_true), "]")

# Vector de validación incluyendo lambda
betas_true_all <- c(betas_true, "lambda" = LAMBDA_TRUE)

# Configuración Stan
params_stan <- list(
  n_chains = N_CORES, 
  n_warmup = 5000,
  n_iter   = 55000,
  refresh  = 0
)

# Compilación del modelo
cat("Compilando modelo Stan...\n")
compiled_model <- rstan::stan_model(file = here("models/power_logit.stan"))

# --- Loop de simulación ---
set.seed(123) 
sim_seeds <- sample.int(1e6, N_SIM)
results_list <- list()

pb <- progress_bar$new(
  format = " Simulando [:bar] :percent | ETA: :eta",
  total = N_SIM, width = 60
)

for (i in 1:N_SIM) {
  
  # Generar datos y estructura para Stan
  data_sim <- simulate_power_logit_data(
    n = N_OBS, betas = betas_true, lambda = LAMBDA_TRUE, 
    rho = RHO, seed = sim_seeds[i], scale_X = TRUE
  )
  
  stan_data <- list(
    N = N_OBS, P = length(betas_true),
    X = data_sim$X, y = as.vector(data_sim$y), sigma = 10
  )
  
  # Ejecución del driver
  res <- run_stan_power_logit_parallel(
    data_list = stan_data, params = params_stan, 
    model_obj = compiled_model, base_seed = sim_seeds[i]
  )
  
  if (!is.null(res$resumen)) {
    df <- res$resumen
    
    # Metadatos y diagnósticos
    df$iteration  <- i
    df$scenario   <- SCENARIO_ID
    df$time_total <- res$times$total
    df$multi_ess  <- res$global_metrics$multivariate_ess
    
    # Métricas de precisión y cobertura
    df$beta_true <- betas_true_all[match(df$parametro, names(betas_true_all))]
    df$bias      <- df$mean - df$beta_true
    df$sq_error  <- df$bias^2
    df$coverage  <- as.integer(df$beta_true >= df$lower_ci & df$beta_true <= df$upper_ci)
    
    results_list[[i]] <- df
    
    # Exportar muestras
    saveRDS(res$muestras, 
            file = sprintf("output/chains/stan/%s_iter_%03d.rds", SCENARIO_ID, i), 
            compress = "xz")
  }
  
  if (i %% 10 == 0) gc()
  pb$tick()
}

# --- Exportación y reporte ---
final_df <- do.call(rbind, results_list)
write.csv(final_df, sprintf("output/%s_results.csv", SCENARIO_ID), row.names = FALSE)

cat("\nProceso finalizado. Resultados en:", SCENARIO_ID, "\n")

# Resumen de recuperación de parámetros
summary_stats <- final_df %>%
  group_by(parametro) %>%
  summarise(
    Bias = mean(bias),
    RMSE = sqrt(mean(sq_error)),
    Cov  = mean(coverage),
    Rhat = mean(rhat),
    ESS  = mean(ess_bulk)
  )

cat("\n--- RESULTADOS: LAMBDA ---\n")
print(summary_stats %>% filter(parametro == "lambda"))

cat("\n--- RESULTADOS: BETAS (Top 5) ---\n")
print(head(summary_stats %>% filter(grepl("beta", parametro)), 5))