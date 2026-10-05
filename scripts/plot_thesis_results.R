# ==============================================================================
# SCRIPT DE VISUALIZACIÓN PARA TESIS (CAPÍTULO APLICACIÓN)
# Genera gráficos de alta calidad a partir de los outputs del pipeline
# ==============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork) # Para combinar gráficos (install.packages("patchwork"))
  library(stringr)
})

# Definir un TEMA ACADÉMICO limpio para todos los gráficos
theme_thesis <- theme_bw() +
  theme(
    text = element_text(family = "sans", color = "black"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold")
  )

# ==============================================================================
# GRÁFICO 1: TRAZA DE CONVERGENCIA (PL vs RPL)
# Réplica mejorada de 'fig_traceplots_comparativos.png'
# ==============================================================================
cat(">>> Generando Gráfico 1: Trazas de Convergencia...\n")

# Cargar datos
df_trace <- readRDS(here("output", "plots", "lambda_trace_data.rds"))

# Procesamiento para el gráfico
# 1. Crear log(Lambda) para visualizar mejor (especialmente RPL que diverge)
# 2. Filtrar iteraciones (thining visual) si son demasiadas para que no pese el PDF
df_plot_trace <- df_trace %>%
  mutate(
    LogLambda = log(Lambda),
    Modelo_Label = ifelse(Modelo == "POWER", 
                          "A) Modelo Power Logit (PL)\nConvergencia Robusta", 
                          "B) Modelo Reversal Power Logit (RPL)\nFallo en Convergencia")
  ) %>%
  filter(Iteracion %% 10 == 0) # Plotear 1 de cada 10 puntos para aligerar

# Colores personalizados (Azules para PL, Rojos para RPL)
p_trace <- ggplot(df_plot_trace, aes(x = Iteracion, y = LogLambda, group = Cadena, color = Cadena)) +
  geom_line(alpha = 0.6, size = 0.3) +
  facet_wrap(~Modelo_Label, scales = "free_y") + # Escalas libres porque RPL se dispara
  scale_color_brewer(palette = "Paired") +
  labs(
    title = "Diagnóstico de Mezcla: Parámetro de Asimetría",
    subtitle = "Comparación de las trazas de log(lambda) entre modelos asimétricos",
    x = "Iteraciones (Post-Warmup)",
    y = expression(log(lambda))
  ) +
  theme_thesis +
  theme(legend.position = "none") # No necesitamos leyenda de cadenas

# Guardar
ggsave(here("output", "plots", "fig1_traza_lambda_comparativa.pdf"), p_trace, width = 10, height = 5)
ggsave(here("output", "plots", "fig1_traza_lambda_comparativa.png"), p_trace, width = 10, height = 5, dpi = 300)


# ==============================================================================
# GRÁFICO 2: FOREST PLOT (COEFICIENTES + LAMBDA)
# Réplica mejorada de 'fig_forest_plot_coeficientes_lambda.png'
# ==============================================================================
cat(">>> Generando Gráfico 2: Forest Plot (Coeficientes)...\n")

# Cargar datos
df_est <- read.table(here("output", "tables", "thesis_power_logit_estimates.csv"), sep = ";", header = TRUE)

# Separar Lambda de Betas
df_lambda <- df_est %>% filter(parametro == "lambda")
df_betas  <- df_est %>% filter(parametro != "lambda")

# Ordenar Betas por magnitud de la media (para que se vea ordenado como tu imagen)
df_betas <- df_betas %>%
  mutate(parametro = reorder(parametro, mean))

# --- Panel A: Betas ---
p_betas <- ggplot(df_betas, aes(x = mean, y = parametro)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.3, color = "navy", size = 0.8) +
  geom_point(size = 3, color = "navy") +
  labs(
    title = "Estimaciones Posteriores: Power Logit",
    subtitle = "A) Coeficientes de Regresión (Escala Log-Odds)",
    x = NULL, # Quitamos X label para unirlo con el de abajo
    y = NULL
  ) +
  theme_thesis +
  theme(panel.grid.major.y = element_line(color = "gray90"))

# --- Panel B: Lambda ---
p_lambda <- ggplot(df_lambda, aes(x = mean, y = parametro)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 0.8) + # Ref = 1 (Logit)
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2, color = "darkgreen", size = 0.8) +
  geom_point(size = 3, color = "darkgreen") +
  labs(
    subtitle = "B) Parámetro de Forma (Escala Natural)",
    x = "Estimación (Media Posterior e IC 95%)",
    y = NULL
  ) +
  xlim(0, max(df_lambda$upper_ci) * 1.1) + # Asegurar que empiece en 0
  theme_thesis +
  theme(
    axis.text.y = element_text(face = "bold.italic", size = 12)
  )

# Combinar con patchwork (Betas ocupa 3/4 del espacio, Lambda 1/4)
p_forest <- p_betas / p_lambda + plot_layout(heights = c(3, 1))

# Guardar
ggsave(here("output", "plots", "fig2_forest_plot_estimates.pdf"), p_forest, width = 9, height = 7)
ggsave(here("output", "plots", "fig2_forest_plot_estimates.png"), p_forest, width = 9, height = 7, dpi = 300)


# ==============================================================================
# GRÁFICO 3: DENSIDAD POSTERIOR DE LAMBDA (CORREGIDO)
# Solución: Usar un dataframe resumen para evitar pintar 398k puntos
# ==============================================================================
cat(">>> Generando Gráfico 3: Densidad Posterior Lambda...\n")

# 1. Preparar datos de las trazas
df_lambda_draws <- df_trace %>% 
  filter(Modelo == "POWER") %>%
  filter(Iteracion > 500)

# 2. Crear un dataframe PEQUEÑO (1 fila) solo para las estadísticas
# Esto evita el warning y hace el gráfico más ligero
stats_lambda <- data.frame(
  media = mean(df_lambda_draws$Lambda),
  li    = quantile(df_lambda_draws$Lambda, 0.025),
  ls    = quantile(df_lambda_draws$Lambda, 0.975),
  y_pos = 0.01 # Altura de la barra
)

# 3. Graficar
p_density <- ggplot(df_lambda_draws, aes(x = Lambda)) +
  # Capa de densidad usa los datos completos (398k filas)
  geom_density(fill = "darkgreen", alpha = 0.3, color = "darkgreen") +
  
  # Línea de referencia Logit
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 1.5, y = 0, label = "Ref: Logit (Simétrico)", 
           color = "red", angle = 90, vjust = -0.5, hjust = 0, size = 3.5) +
  
  # CORRECCIÓN: Usar 'data = stats_lambda' y 'inherit.aes = FALSE'
  # Así solo dibuja 1 barra y 1 punto, en lugar de 398,000.
  geom_errorbarh(data = stats_lambda, 
                 aes(xmin = li, xmax = ls, y = y_pos), 
                 height = 0.02, color = "black", size = 1, inherit.aes = FALSE) +
  
  geom_point(data = stats_lambda, 
             aes(x = media, y = y_pos), 
             size = 3, inherit.aes = FALSE) +
  
  labs(
    title = "Evidencia de Asimetría en el Fraude",
    subtitle = "Distribución Posterior del parámetro Lambda (Power Logit)",
    x = expression(lambda ~ "(Parámetro de Forma)"),
    y = "Densidad"
  ) +
  theme_thesis

# Guardar (Ajusta la ruta "plots" o "figures" según tu carpeta creada)
ggsave(here("output", "plots", "fig3_densidad_lambda.pdf"), p_density, width = 8, height = 5)
