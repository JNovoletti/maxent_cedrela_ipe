# =========================================================================
# Objetivo: Extrair variáveis ambientais de rasters (.tif) para os pontos
#           de ocorrência das espécies Cedrela e Handroanthus.
# Baseado em: 02_extr_dados_raster.qmd (Prof. Dr. Deoclecio Jardim Amorim)
# =========================================================================

# ---------------------------------------------------------------------------
# 1. Configuração inicial
# ---------------------------------------------------------------------------

# Limpar ambiente
# rm(list = ls())
# gc(reset = TRUE)
# graphics.off()

# Pacotes necessários
if (!require(readxl))    install.packages("readxl",    dep = TRUE, quiet = TRUE)
if (!require(writexl))   install.packages("writexl",   dep = TRUE, quiet = TRUE)
if (!require(tidyverse)) install.packages("tidyverse", dep = TRUE, quiet = TRUE)
if (!require(terra))     install.packages("terra",     dep = TRUE, quiet = TRUE)
if (!require(sf))        install.packages("sf",        dep = TRUE, quiet = TRUE)
if (!require(here))      install.packages("here",      dep = TRUE, quiet = TRUE)

library(readxl)
library(writexl)
library(tidyverse)
library(terra)
library(sf)
library(here)

# ---------------------------------------------------------------------------
# 2. Caminhos
# ---------------------------------------------------------------------------

# Pasta dos rasters (caminho absoluto conforme solicitado)
raster_dir <- "C:/Users/Aluno/Documents/dslab/Rasters"

# Pasta de dados de entrada/saída (relativa ao projeto)
data_dir <- here("data")

# Garantir que a pasta de saída exista
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# 3. Carregamento e empilhamento dos rasters
# ---------------------------------------------------------------------------

# Listar arquivos .tif na pasta, excluindo o stack já consolidado
raster_files <- list.files(
  path       = raster_dir,
  pattern    = "\\.tif$",
  full.names = TRUE
)

# Excluir arquivos que não são variáveis ambientais individuais:
excluir_pattern <- "Brazil_env_stack|Amazon_stack|Caatinga_elev_corridor|elev_corridor_minus1"
raster_files <- raster_files[!grepl(excluir_pattern, raster_files)]

cat("Total de rasters encontrados:", length(raster_files), "\n")

# Empilhar como SpatRaster multicamadas
r_stack <- terra::rast(raster_files)

# Nomear camadas a partir dos nomes dos arquivos (sem extensão)
names(r_stack) <- tools::file_path_sans_ext(basename(raster_files))

cat("Camadas no stack:", nlyr(r_stack), "\n")
cat("CRS do stack:", terra::crs(r_stack, describe = TRUE)$name, "\n")

# ---------------------------------------------------------------------------
# 4. Função genérica para extrair variáveis para uma espécie
# ---------------------------------------------------------------------------

extrair_variaveis <- function(csv_path, especie, r_stack, out_dir) {
  
  cat("\n--- Processando:", especie, "---\n")
  
  # 4.1 Leitura dos dados
  dados <- read.csv2(csv_path, header = TRUE, stringsAsFactors = FALSE)
  cat("Linhas lidas:", nrow(dados), "\n")
  cat("Colunas disponíveis:", paste(colnames(dados), collapse = ", "), "\n")
  
  # 4.2 Identificar colunas de coordenadas de forma robusta
  #     (aceita latitude/longitude, lat/lon, lat/long, decimalLatitude/decimalLongitude)
  col_lon <- intersect(
    tolower(colnames(dados)),
    c("longitude", "lon", "long", "decimallongitude", "x")
  )[1]
  
  col_lat <- intersect(
    tolower(colnames(dados)),
    c("latitude", "lat", "decimallatitude", "y")
  )[1]
  
  if (is.na(col_lon) || is.na(col_lat)) {
    stop("Não foi possível identificar as colunas de latitude/longitude em ", csv_path)
  }
  
  # Recuperar nomes originais (case-sensitive)
  col_lon_orig <- colnames(dados)[tolower(colnames(dados)) == col_lon]
  col_lat_orig <- colnames(dados)[tolower(colnames(dados)) == col_lat]
  
  par(mfrow = c(1,1))
  plot(resultados$lon, resultados$lat, col = "darkgreen", pch = 19, cex = 0.5, xlab = "Longitude", ylab = "Latitude")
  maps::map(add = TRUE)
  
  cat("Coluna de longitude:", col_lon_orig, "\n")
  cat("Coluna de latitude:",  col_lat_orig, "\n")
  
  # 4.3 Remover linhas sem coordenadas válidas
  dados <- dados %>%
    dplyr::filter(
      !is.na(.data[[col_lon_orig]]),
      !is.na(.data[[col_lat_orig]])
    )
  cat("Linhas com coordenadas válidas:", nrow(dados), "\n")
  
  # 4.4 Converter para objeto sf (WGS84 EPSG:4674 - SIRGAS 2000)
  dados_sf <- sf::st_as_sf(
    dados,
    coords = c(col_lon_orig, col_lat_orig),
    crs    = 4674,
    remove = FALSE
  )
  
  # 4.5 Reprojetar para o CRS do raster, se necessário
  if (!sf::st_crs(dados_sf) == sf::st_crs(r_stack)) {
    cat("Reprojetando pontos para o CRS dos rasters...\n")
    dados_sf <- sf::st_transform(dados_sf, sf::st_crs(r_stack))
  }
  
  cat("CRS dos rasters:\n")
  print(terra::crs(r_stack, describe = TRUE))
  
  cat("\nCRS dos pontos:\n")
  print(sf::st_crs(dados_sf, describe = TRUE))
  
  # 4.6 Extração dos valores (método bilinear: interpola entre as 4 células vizinhas)
  cat("Extraindo valores raster (pode demorar)...\n")
  extracted_values <- terra::extract(r_stack, terra::vect(dados_sf), method = "bilinear")
  
  # Remover a coluna ID gerada pela função extract
  extracted_values <- extracted_values[, -1, drop = FALSE]
  
  # 4.7 Combinar dados originais + valores extraídos
  coords <- sf::st_coordinates(dados_sf)
  resultado <- dados_sf %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(x = coords[, 1], y = coords[, 2]) %>%
    dplyr::relocate(x, y) %>%
    dplyr::bind_cols(extracted_values)
  
  # 4.8 Exportação
  out_xlsx <- file.path(out_dir, paste0(especie, "_var_amb.xlsx"))
  out_csv  <- file.path(out_dir, paste0(especie, "_var_amb.csv"))
  
  writexl::write_xlsx(resultado, out_xlsx)
  write.csv(resultado, out_csv, row.names = FALSE)
  
  cat("Arquivo salvo em:", out_xlsx, "\n")
  cat("Arquivo salvo em:", out_csv,  "\n")
  cat("Dimensões finais:", nrow(resultado), "x", ncol(resultado), "\n")
  
  return(resultado)
}

# ---------------------------------------------------------------------------
# 5. Execução para cada espécie
# ---------------------------------------------------------------------------

cedrela_result <- extrair_variaveis(
  csv_path = file.path(data_dir, "cedrela_br.csv"),
  especie  = "cedrela_br",
  r_stack  = r_stack,
  out_dir  = data_dir
)

handroanthus_result <- extrair_variaveis(
  csv_path = file.path(data_dir, "handroanthus_br.csv"),
  especie  = "handroanthus",
  r_stack  = r_stack,
  out_dir  = data_dir
)

print(head(cedrela_result))

print(head(handroanthus_result))
