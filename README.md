# MCMC Thesis: MALA, Stan y JAGS

Este repositorio contiene el entorno de trabajo y los scripts de simulación desarrollados para mi tesis de maestría. Se utiliza Docker para garantizar la reproducibilidad de los resultados, gestionando las dependencias de JAGS y los compiladores de C++ necesarios para RTMB y Stan.

## Resumen

Tesis de la Maestría en Estadística de la PUCP (2026), aprobada con calificación sobresaliente. Título: *Fundamentos y aplicaciones de las ecuaciones de Langevin en métodos MCMC*. Asesor: Dr. Luis Valdivieso.

- **Teoría:** difusión de Langevin, discretización de Euler–Maruyama (ULA), el algoritmo MALA (*Metropolis-Adjusted Langevin Algorithm*) y la ergodicidad mediante cadenas de Harris.
- **Implementación:** MALA desde cero en R, con diferenciación automática (RTMB) para el gradiente y ajuste estocástico del tamaño de paso (Robbins–Monro).
- **Estudio de simulación:** comparación con JAGS y Stan (NUTS) en regresión logística simétrica (Logit) y asimétrica (Power Logit), con distintos niveles de dimensión y correlación.
  - En escenarios convexos estándar, MALA escala mejor que JAGS.
  - Con alta dimensión y fuerte correlación o asimetría, necesita muchas más iteraciones que Stan para lograr un tamaño de muestra efectivo competitivo.
- **Aplicación:** detección de fraude en subastas en línea (*shill bidding*), un evento raro. El enlace asimétrico (Power Logit) supera claramente al logit estándar.

## Requisitos
- Docker Desktop
- Git

## Configuración del Entorno

1. Clonar el repositorio

```
git clone git@github.com:bxcalsisilva/msc-thesis-mala.git
cd msc-thesis-mala
```

2. Construcción de la imagen

La imagen utiliza R e incluye las herramientas de compilación optimizadas para modelos asimétricos.

```
docker build -t msc-thesis-env .
```

3. Ejecución

El script de PowerShell adjunto automatiza la limpieza de sesiones, el montaje de volúmenes y la apertura de RStudio.

```
./start_thesis.ps1
```

## Acceso a RStudio
El entorno estará disponible en http://localhost:8787 una vez el contenedor esté activo.

- Usuario: `rstudio`
- Contraseña: `thesis`

## Estructura del Proyecto
- `R/`: Funciones core y drivers para los muestreadores. 
- `models/`: Definición de modelos en código Stan y JAGS. 
- `scripts/`: Scripts de ejecución para simulaciones y aplicación. 
- `data/`: Dataset de Shill Bidding utilizado en la tesis. 
- `output/`: Tablas de resultados, métricas de convergencia y gráficos. 
