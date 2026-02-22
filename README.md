# MCMC Thesis: MALA, Stan y JAGS

Este repositorio contiene el entorno de trabajo y los scripts de simulación desarrollados para mi tesis de maestría. Se utiliza Docker para garantizar la reproducibilidad de los resultados, gestionando las dependencias de JAGS y los compiladores de C++ necesarios para RTMB y Stan.

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
