# 1. Detener el contenedor si ya estaba corriendo (para evitar errores)
# Write-Host "Clearning previous session..."
docker stop contenedor_tesis 2>$null

# 2. Iniciar el contenedor (El comando largo)
# Write-Host "Starting thesis environment"
docker run --rm -d `
  -p 8787:8787 `
  -e PASSWORD=thesis `
  -e ROOT=true `
  -v ${PWD}:/home/rstudio/msc-thesis-mala `
  --name contenedor_tesis `
  msc-thesis-env

# 3. Esperar unos segundos a que RStudio arranque
# Write-Host "Waiting for Rstudio"
Start-Sleep -Seconds 3

# 4. Abrir el navegador automáticamente
# Write-Host "Starting Rstudio"
Start-Process "http://localhost:8787"