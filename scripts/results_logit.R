# scripts/process_results.R
rm(list = ls())
library(tidyverse)
library(here)

# --- Configuración de rutas ---
path_files <- here("output")

if (!dir.exists(path_files)) {
  stop("Carpeta 'output' no encontrada.")
}

files <- list.files(path = path_files, 
                    pattern = "_results\\.csv$", 
                    full.names = TRUE)

# --- Función de procesamiento ---
read_and_process <- function(filepath) {
  fname <- basename(filepath)
  
  # Extracción de metadatos desde el nombre del archivo
  algo       <- str_extract(fname, "^[a-z]+") 
  is_power   <- str_detect(fname, "power")
  model_type <- ifelse(is_power, "Power Logit", "Logit Estándar")
  
  n_val   <- as.numeric(str_extract(fname, "(?<=_n)[0-9]+"))
  p_val   <- as.numeric(str_extract(fname, "(?<=_p)[0-9]+"))
  rho_val <- as.numeric(str_extract(fname, "(?<=_rho)[0-9]+")) / 10
  
  lambda_val <- NA
  if (is_power) {
    lambda_val <- as.numeric(str_extract(fname, "(?<=_s)[0-9]+")) / 10
  }
  
  # Carga y estructuración
  read_csv(filepath, show_col_types = FALSE) %>%
    mutate(
      algorithm   = toupper(algo),
      model       = model_type,
      n           = n_val,
      p           = p_val,
      rho         = rho_val,
      lambda_sim  = lambda_val,
      source_file = fname
    ) %>%
    relocate(algorithm, model, n, p, rho, lambda_sim)
}

# --- Ejecución de la lectura ---
if (length(files) > 0) {
  cat("Procesando", length(files), "archivos de resultados...\n")
  master_df <- map_dfr(files, read_and_process)
} else {
  stop("No se encontraron archivos '_results.csv' en la carpeta output.")
}

# --- Análisis: Modelo Logit Estándar ---
logit_df <- master_df %>% filter(model == 'Logit Estándar')

# Agregación por réplica
replica_stats <- logit_df %>%
  group_by(algorithm, n, p, rho, iteration) %>%
  summarise(
    max_rhat_rep   = max(rhat, na.rm = TRUE),
    avg_mcse_ratio = mean(mcse_mean / sd, na.rm = TRUE),
    .groups = "drop"
  )

# Tabla resumen de convergencia
logit_conv_table <- replica_stats %>%
  group_by(algorithm, n, p, rho) %>%
  summarise(
    Avg_Max_Rhat = mean(max_rhat_rep, na.rm = TRUE),
    Worst_Rhat   = max(max_rhat_rep, na.rm = TRUE),
    Conv_Rate    = mean(max_rhat_rep < 1.05, na.rm = TRUE) * 100,
    Avg_MCSE     = mean(avg_mcse_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(n, p, rho, algorithm)

# --- Exportación ---
print(logit_conv_table)

output_path <- here("output", "tables", "logit_convergence_summary.csv")
if(!dir.exists(dirname(output_path))) dir.create(dirname(output_path), recursive = TRUE)

write_csv(logit_conv_table, output_path)
cat("\nResumen guardado en:", output_path, "\n")