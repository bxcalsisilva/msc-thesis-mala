# 1. Base estable
FROM rocker/rstudio:4.4.2

# 2. Instalacion de dependencias del sistema
RUN apt-get update && apt-get install -y --no-install-recommends \
    jags \
    build-essential \
    libfftw3-dev \
    libglpk-dev \
    libxml2-dev \
    libcairo2-dev \
    libgit2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# 3. Configuracion de Compiladores (Makevars)
# Esto le dice a R: "Usa g++ con el estandar C++14 y C++17, y optimiza el codigo al maximo (-O3)"
# Usamos un solo RUN para crear el directorio y el archivo de configuracion
RUN mkdir -p /home/rstudio/.R && \
    echo "CXX14FLAGS=-O3 -march=native -mtune=native -fPIC" >> /home/rstudio/.R/Makevars && \
    echo "CXX14=g++" >> /home/rstudio/.R/Makevars && \
    echo "CXX17FLAGS=-O3 -march=native -mtune=native -fPIC" >> /home/rstudio/.R/Makevars && \
    echo "CXX17=g++" >> /home/rstudio/.R/Makevars && \
    chown -R rstudio:rstudio /home/rstudio/.R

# 4. Instalacion de paquetes de R (Ordenados por importancia)
# Nota: Agrupamos rstan y RTMB primero porque son los mas pesados
RUN install2.r --error --skipinstalled \
    rjags \
    runjags \
    rstan \
    RTMB \
    tidyverse \
    dplyr \
    here \
    progress \
    mcmcse \
    stringr 

# CORRER: docker build -t msc-thesis-env .