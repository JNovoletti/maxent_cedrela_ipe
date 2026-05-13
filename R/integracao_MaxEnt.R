# =========================================================================
# Objetivo: Avaliar modelos de nicho ecológico para Cedrela e Handroanthus 
#           no Brasil utilizando o algoritmo MaxEnt
# Baseado em: "Modelling the Effects of Climate Change Using Maxent and R"
#             (Hidasi-Neto)
# =========================================================================

#'
## -----------------------------------------------------------------------------
#| include: false

# Ajuste global para gráficos
library(knitr)
opts_knit$set(global.par = TRUE)
par(mar = c(5, 5, 1, 1))

# ---------------------------------------------------------------------------
# 1. Configuração inicial
# ---------------------------------------------------------------------------

# Limpar ambiente
# rm(list = ls())
# gc(reset = TRUE)
# graphics.off()

# Pacotes necessários
if(!require(readxl))      install.packages("readxl",      dep = TRUE, quiet = TRUE)
if(!require(writexl))     install.packages("writexl",     dep = TRUE, quiet = TRUE)
if(!require(tidyverse))   install.packages("tidyverse",   dep = TRUE, quiet = TRUE)
if(!require(terra))       install.packages("terra",       dep = TRUE, quiet = TRUE)
if(!require(sf))          install.packages("sf",          dep = TRUE, quiet = TRUE)
if(!require(sp))          install.packages("sp",          dep = TRUE, quiet = TRUE)
if(!require(here))        install.packages("here",        dep = TRUE, quiet = TRUE)
if(!require(dismo))       install.packages("dismo",       dep = TRUE, quiet = TRUE)
if(!require(raster))      install.packages("raster",      dep = TRUE, quiet = TRUE)
if(!require(rJava))       install.packages("rJava",       dep = TRUE, quiet = TRUE)
if(!require(maps))        install.packages("maps",        dep = TRUE, quiet = TRUE)
if(!require(geobr))       install.packages("geobr",       dep = TRUE, quiet = TRUE)

# Observação importante:
# O pacote `dismo` requer o arquivo `maxent.jar` (binário Java).
# Caso não venha junto:
# 1. Baixe em: https://biodiversityinformatics.amnh.org/open_source/maxent/
# 2. Copie `maxent.jar` para: system.file("java", package = "dismo")
# 3. Confirme com: file.exists(paste0(system.file("java", package = "dismo"),
#                                     "/maxent.jar"))
maxent_jar <- file.path(system.file("java", package = "dismo"), "maxent.jar")
if (!file.exists(maxent_jar)) {
  warning("maxent.jar nao encontrado em: ", maxent_jar,
          "\nBaixe em https://biodiversityinformatics.amnh.org/open_source/maxent/")
}


#'
#' # 2. Reprodutibilidade
#'
## ----seed---------------------------------------------------------------------
set_reproducibility <- function(seed = 1350) {
  set.seed(seed)
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")
}

# Fixando a reprodutibilidade
set_reproducibility()


#'
#' # 3. Carregando os conjuntos de dados
#'
## ----importacao_dados---------------------------------------------------------
# Dados de Cedrela
cedrela <- readxl::read_excel(here::here("data", "cedrela_br_var_amb.xlsx"), sheet = 1)
head(cedrela)

# Dados de Handroanthus
handroanthus <- readxl::read_excel(here::here("data", "handroanthus_var_amb.xlsx"), sheet = 1)
head(handroanthus)


#'
#' # 4. Carregamento do stack de rasters
#'
#' MaxEnt precisa do **raster multicamadas completo** para gerar o background
#' (pseudo-ausências) e para predizer a adequabilidade em toda a área de estudo.
#'
## ----stack--------------------------------------------------------------------
# Caminho dos rasters
raster_dir <- "C:/Users/Aluno/Documents/dslab/Rasters"

# Listar arquivos .tif
# raster_files <- list.files(path = raster_dir, pattern = "\\.tif$", full.names = TRUE)
raster_files <- list.files(  path = raster_dir,   pattern = "32\\.vapor|33\\.solar\\.tif$",   full.names = TRUE)

# Excluir stacks consolidados e recortes regionais
excluir_pattern <- "Brazil_env_stack|Amazon_stack|Caatinga_elev_corridor|elev_corridor_minus1"
raster_files <- raster_files[!grepl(excluir_pattern, raster_files)]

# Empilhar como SpatRaster
r_stack <- terra::rast(raster_files)
names(r_stack) <- tools::file_path_sans_ext(basename(raster_files))

cat("Camadas no stack:", nlyr(r_stack), "\n")
cat("Extensao espacial:\n"); print(terra::ext(r_stack))

# Para usar com dismo/maxent, precisamos do objeto RasterStack do pacote raster
r_stack_raster <- raster::stack(raster_files)
names(r_stack_raster) <- tools::file_path_sans_ext(basename(raster_files))


#'
#' # 5. Preparação dos dados de ocorrência
#'
#' Para cada espécie, vamos:
#'
#' - Extrair apenas as coordenadas (longitude, latitude);
#' - Remover duplicatas e registros sem coordenada;
#' - Dividir em **treino (75%)** e **teste (25%)** usando k-fold.
#'
## ----preparacao---------------------------------------------------------------
preparar_ocorrencias <- function(dados, especie) {
  
  cat("\n--- Preparando:", especie, "---\n")
  
  # Padronizar nomes das colunas de coordenadas (script 02 cria x = lon, y = lat)
  if (all(c("x", "y") %in% colnames(dados))) {
    pontos <- data.frame(lon = dados$x, lat = dados$y)
  } else if (all(c("longitude", "latitude") %in% colnames(dados))) {
    pontos <- data.frame(lon = dados$longitude, lat = dados$latitude)
  } else {
    stop("Colunas de coordenadas nao encontradas para ", especie)
  }
  
  # Remover NAs e duplicatas
  pontos <- pontos[complete.cases(pontos), ]
  pontos <- pontos[!duplicated(pontos), ]
  
  cat("Pontos validos:", nrow(pontos), "\n")
  
  return(pontos)
}

# Aplicar para as duas especies
cedrela_pts      <- preparar_ocorrencias(cedrela,      "Cedrela")
handroanthus_pts <- preparar_ocorrencias(handroanthus, "Handroanthus")

# Visualizar distribuicao geografica
par(mfrow = c(1, 2))
plot(cedrela_pts$lon, cedrela_pts$lat,
     col = "darkgreen", pch = 19, cex = 0.5,
     xlab = "Longitude", ylab = "Latitude", main = "Cedrela")
maps::map(add = TRUE)

plot(handroanthus_pts$lon, handroanthus_pts$lat,
     col = "darkorange", pch = 19, cex = 0.5,
     xlab = "Longitude", ylab = "Latitude", main = "Handroanthus")
maps::map(add = TRUE)
par(mfrow = c(1, 1))


#'
#' # 6. Funcao para ajuste do modelo MaxEnt
#'
#' Encapsulamos o pipeline MaxEnt em uma funcao reutilizavel para garantir que
#' as duas especies passem exatamente pelo mesmo procedimento.
#'
## ----funcao_maxent------------------------------------------------------------
ajustar_maxent <- function(pontos, climate_stack, especie,
                           ext, n_background = 1000, k = 5) {
  
  cat("\n========== MaxEnt:", especie, "==========\n")
  
  # 6.1 Divisao treino (80%) / teste (20%) com k-fold
  set_reproducibility()
  group      <- dismo::kfold(pontos, k)
  pres_train <- pontos[group != 1, ]
  pres_test  <- pontos[group == 1, ]
  
  cat("Treino:", nrow(pres_train), "| Teste:", nrow(pres_test), "\n")
  
  # 6.2 Ajuste MaxEnt (requer maxent.jar)
  cat("Ajustando MaxEnt...\n")
  xm <- dismo::maxent(climate_stack, pres_train)
  
  # 6.3 Geracao de background (pseudo-ausencias)
  set_reproducibility()
  backg <- dismo::randomPoints(climate_stack, n = n_background,
                               ext = ext, extf = 1.25)
  colnames(backg) <- c("lon", "lat")
  
  group_bg   <- dismo::kfold(backg, k)
  backg_train <- backg[group_bg != 1, ]
  backg_test  <- backg[group_bg == 1, ]
  
  # 6.4 Avaliacao do modelo (AUC)
  e <- dismo::evaluate(pres_test, backg_test, xm, climate_stack)
  cat("AUC:", round(e@auc, 4), "\n")
  
  # 6.5 Predicao espacial (mapa de adequabilidade)
  cat("Gerando predicao espacial...\n")
  p_pred <- dismo::predict(climate_stack, xm, ext = ext, progress = "")
  
  # Retornar lista com todos os objetos relevantes
  return(list(
    especie     = especie,
    modelo      = xm,
    avaliacao   = e,
    auc         = e@auc,
    pres_train  = pres_train,
    pres_test   = pres_test,
    backg_train = backg_train,
    backg_test  = backg_test,
    predicao    = p_pred
  ))
}


#'
#' # 7. Ajuste dos modelos para cada especie
#'
## ----ajuste-cedrela-----------------------------------------------------------
resultado_cedrela <- ajustar_maxent(
  pontos        = cedrela_pts,
  climate_stack = r_stack_raster,
  especie       = "Cedrela",
  ext           = ext_brasil,
  n_background  = 1000,
  k             = 5
)

#'
## ----ajuste-handroanthus------------------------------------------------------
resultado_handroanthus <- ajustar_maxent(
  pontos        = handroanthus_pts,
  climate_stack = r_stack_raster,
  especie       = "Handroanthus",
  ext           = ext_brasil,
  n_background  = 1000,
  k             = 5
)


#'
#' # 8. Importancia das variaveis ambientais
#'
#' MaxEnt fornece a contribuicao percentual de cada variavel para o modelo.
#' Isso e analogo ao Variable Importance Plot (VIP) do Random Forest.
#'
## ----importancia--------------------------------------------------------------
# Cedrela
par(mfrow = c(1, 2))
plot(resultado_cedrela$modelo,
     main = "Importancia das variaveis - Cedrela")

# Handroanthus
plot(resultado_handroanthus$modelo,
     main = "Importancia das variaveis - Handroanthus")
par(mfrow = c(1, 1))


#'
#' # 9. Avaliacao dos modelos - Curvas ROC
#'
#' A curva ROC mostra o trade-off entre sensibilidade e especificidade.
#' AUC > 0.7 e considerado bom; AUC > 0.9 e considerado excelente.
#'
## ----roc----------------------------------------------------------------------
par(mfrow = c(1, 2))
plot(resultado_cedrela$avaliacao, "ROC",
     main = paste0("Cedrela (AUC = ", round(resultado_cedrela$auc, 3), ")"))

plot(resultado_handroanthus$avaliacao, "ROC",
     main = paste0("Handroanthus (AUC = ", round(resultado_handroanthus$auc, 3), ")"))
par(mfrow = c(1, 1))


#'
#' # 10. Mapas de adequabilidade de habitat (presente)
#'
#' Visualizamos os mapas preditivos para as duas especies, com pontos de
#' ocorrencia sobrepostos.
#'
## ----mapas-presente-----------------------------------------------------------
# Carregar contorno do Brasil
brasil <- geobr::read_country(year = 2020)

# Mapa Cedrela
par(mfrow = c(1, 2))
plot(resultado_cedrela$predicao,
     main = "Adequabilidade - Cedrela",
     col = hcl.colors(100, "YlGnBu", rev = TRUE))
plot(sf::st_geometry(brasil), add = TRUE, border = "gray30")
points(cedrela_pts$lon, cedrela_pts$lat, pch = 20, cex = 0.4, col = "red")

# Mapa Handroanthus
plot(resultado_handroanthus$predicao,
     main = "Adequabilidade - Handroanthus",
     col = hcl.colors(100, "YlOrRd", rev = TRUE))
plot(sf::st_geometry(brasil), add = TRUE, border = "gray30")
points(handroanthus_pts$lon, handroanthus_pts$lat, pch = 20, cex = 0.4, col = "blue")
par(mfrow = c(1, 1))


#'
#' # 11. Exportacao dos resultados
#'
#' Salvamos os rasters de adequabilidade e os modelos ajustados para uso futuro
#' (por exemplo, projecoes em cenarios de mudanca climatica).
#'
## ----exportacao---------------------------------------------------------------
# Criar pasta de saida
out_dir <- here::here("maxent_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Exportar rasters de predicao
raster::writeRaster(resultado_cedrela$predicao,
                    filename = file.path(out_dir, "maxent_cedrela_presente.tif"),
                    overwrite = TRUE)

raster::writeRaster(resultado_handroanthus$predicao,
                    filename = file.path(out_dir, "maxent_handroanthus_presente.tif"),
                    overwrite = TRUE)

# Salvar modelos R (para reabrir sem refazer o ajuste)
saveRDS(resultado_cedrela,      file.path(out_dir, "modelo_cedrela.rds"))
saveRDS(resultado_handroanthus, file.path(out_dir, "modelo_handroanthus.rds"))

# Tabela-resumo dos AUC
resumo <- data.frame(
  Especie = c("Cedrela", "Handroanthus"),
  AUC     = c(resultado_cedrela$auc, resultado_handroanthus$auc),
  N_treino = c(nrow(resultado_cedrela$pres_train),
               nrow(resultado_handroanthus$pres_train)),
  N_teste  = c(nrow(resultado_cedrela$pres_test),
               nrow(resultado_handroanthus$pres_test))
)
print(resumo)

writexl::write_xlsx(resumo, file.path(out_dir, "resumo_maxent.xlsx"))

