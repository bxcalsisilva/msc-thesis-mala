# scripts/plot_fig_forest.R
# Figura: Forest plot de las estimaciones del Power Logit (MALA).
#   Panel A: coeficientes de regresion (con IC 95%), linea de referencia en 0.
#   Panel B: parametro de forma lambda (escala natural), linea roja en 1.
#
# Lee:  output/tables/mala_shill_summary_power.csv
# Crea: output/plots/fig2_forest_plot_estimates.pdf
#
# NOTA: con el diseno final (4 covariables) el panel A muestra el intercepto
# y las 4 covariables retenidas; con lambda ~ Uniforme(-2,2), el IC de lambda
# queda acotado (aprox. [3.7, 7.4]), a diferencia de la figura del prior normal.
# ==============================================================================
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

ARCHIVO <- here("output", "tables", "mala_shill_summary_power.csv")
SALIDA  <- here("output", "plots", "fig2_forest_plot_estimates.pdf")
dir.create(dirname(SALIDA), showWarnings = FALSE, recursive = TRUE)

COL_PUNTO  <- "#1F3B73"   # azul oscuro (coeficientes)
COL_LAMBDA <- "#1B7837"   # verde (lambda)

# ------------------------------------------------------------------------------
# 1. Datos
# ------------------------------------------------------------------------------
est <- read.csv(ARCHIVO, stringsAsFactors = FALSE)

betas  <- est %>% filter(parametro != "lambda")
lambda <- est %>% filter(parametro == "lambda")

# Ordenar coeficientes de mayor a menor media (el mayor arriba)
betas <- betas %>% mutate(parametro = factor(parametro, levels = parametro[order(mean)]))

# ------------------------------------------------------------------------------
# 2. Panel A: coeficientes de regresion
# ------------------------------------------------------------------------------
pA <- ggplot(betas, aes(x = mean, y = parametro)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.25,
                 color = COL_PUNTO, linewidth = 0.7) +
  geom_point(size = 2.4, color = COL_PUNTO) +
  labs(subtitle = "A) Coeficientes de Regresión (Escala del Predictor Lineal)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(color = COL_PUNTO, size = 10))

# ------------------------------------------------------------------------------
# 3. Panel B: parametro de forma lambda
# ------------------------------------------------------------------------------
pB <- ggplot(lambda, aes(x = mean, y = parametro)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.18,
                 color = COL_LAMBDA, linewidth = 0.8) +
  geom_point(size = 3, color = COL_LAMBDA) +
  scale_y_discrete(labels = expression(lambda)) +
  labs(subtitle = "B) Parámetro de Forma (Escala Natural)",
       x = "Estimación (Media Posterior e IC 95%)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(color = COL_LAMBDA, size = 10))

# ------------------------------------------------------------------------------
# 4. Combinar y guardar
# ------------------------------------------------------------------------------
fig <- (pA / pB) +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(
    title = "Estimaciones Posteriores: Power Logit",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
  )

ggsave(SALIDA, fig, width = 7, height = 5.5, device = cairo_pdf)
cat(sprintf("-> %s\n", SALIDA))
