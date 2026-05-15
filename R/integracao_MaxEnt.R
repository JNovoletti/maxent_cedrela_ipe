# =========================================================================
# Objetivo: Avaliar modelos de nicho ecologico para Cedrela e Handroanthus
#           no Brasil utilizando o algoritmo MaxEnt
# Baseado em: "Modelling the Effects of Climate Change Using Maxent and R"
#             (Hidasi-Neto)
# Selecao de variaveis (VIF + Jackknife) baseada em:
#             Tesfamariam et al. (2022). MaxEnt-based modeling of suitable
#             habitat for rehabilitation of Podocarpus forest at landscape-scale
# =========================================================================

#'
## -----------------------------------------------------------------------------
#| include: false

# Ajuste global para graficos
library(knitr)
opts_knit$set(global.par = TRUE)
par(mar = c(5, 5, 1, 1))

# ---------------------------------------------------------------------------
# 1. Configuracao inicial
# ---------------------------------------------------------------------------

# Limpar ambiente
# rm(list = ls())
# gc(reset = TRUE)
# graphics.off()

# Pacotes necessarios
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
if(!require(usdm))        install.packages("usdm",        dep = TRUE, quiet = TRUE)
# 'usdm' fornece vif()/vifstep() para diagnostico de colinearidade

# Observacao importante:
# O pacote `dismo` requer o arquivo `maxent.jar` (binario Java).
# Caso nao venha junto:
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

# Fixando a reprodutibilidade
set.seed(1350)


#'
#' # 3. Carregando os conjuntos de dados
#'
## ----importacao_dados---------------------------------------------------------
# Dados de Cedrela
cedrela <- read.csv(here::here("data", "cedrela_br_var_amb.csv"))
head(cedrela)

# Dados de Handroanthus
handroanthus <- read.csv(here::here("data", "handroanthus_var_amb.csv"))
head(handroanthus)


#'
#' # 4. Carregamento do stack de rasters
#'
#' MaxEnt precisa do **raster multicamadas completo** para gerar o background
#' (pseudo-ausencias) e para predizer a adequabilidade em toda a area de estudo.
#'
## ----stack--------------------------------------------------------------------
# Caminho dos rasters

# Listar arquivos .tif
raster_files <- list.files(here::here("rasters"), pattern = "\\.tif$", full.names = TRUE)
#raster_files <- list.files(  path = raster_dir,   pattern = "32\\.vapor|33\\.solar\\.tif$",   full.names = TRUE)

# Excluir stacks consolidados e recortes regionais
excluir_pattern <- "Brazil_env_stack|Amazon_stack|Caatinga_elev_corridor|elev_corridor_minus1|48.density"
raster_files <- raster_files[!grepl(excluir_pattern, raster_files)]

# Empilhar como SpatRaster
r_stack <- terra::rast(raster_files)
names(r_stack) <- tools::file_path_sans_ext(basename(raster_files))

cat("Camadas no stack:", nlyr(r_stack), "\n")
cat("Extensao espacial:\n"); print(terra::ext(r_stack))

# Para usar com dismo/maxent, precisamos do objeto RasterStack do pacote raster
r_stack_raster <- raster::stack(raster_files)
names(r_stack_raster) <- tools::file_path_sans_ext(basename(raster_files))

#Definir a extensao
ext_brasil <- raster::extent(r_stack_raster)


#'
#' # 5. Preparacao dos dados de ocorrencia
#'
#' Para cada especie, vamos:
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
#' # 5.1 Selecao de variaveis ambientais (VIF + Jackknife)
#'
#' Seguindo a metodologia de Tesfamariam et al. (2022), reduzimos o numero de
#' variaveis em duas etapas:
#'   (a) Teste de multicolinearidade via VIF (Variance Inflation Factor),
#'       removendo variaveis com VIF > 5.
#'   (b) Teste de importancia via Jackknife do MaxEnt, removendo variaveis
#'       cuja omissao nao reduz o ganho do modelo (baixa contribuicao).
#'
## ----extrair-valores----------------------------------------------------------
#'
#' ## 5.1.1 Extracao dos valores ambientais nos pontos
#'
#' Para o VIF, precisamos dos valores das variaveis nos pontos de presenca +
#' pontos de background (a colinearidade pode mudar entre regioes do espaco
#' ambiental, por isso incluimos ambos).
#'
extrair_valores_amb <- function(pontos, climate_stack, n_background = 10000,
                                ext, especie) {
  
  cat("\n--- Extraindo valores ambientais:", especie, "---\n")
  
  # Gerar background com a mesma extensao usada no MaxEnt
  set.seed(1350)
  backg <- dismo::randomPoints(climate_stack, n = n_background,
                               ext = ext, extf = 1.25)
  colnames(backg) <- c("lon", "lat")
  
  # Extrair valores das variaveis ambientais
  vals_pres <- raster::extract(climate_stack, pontos, method="bilinear")
  vals_back <- raster::extract(climate_stack, backg, method="bilinear")
  
  # Unir presencas + background para o diagnostico
  vals_all <- rbind(vals_pres, vals_back)
  vals_all <- as.data.frame(vals_all)
  vals_all <- vals_all[complete.cases(vals_all), ]
  
  cat("Linhas validas para diagnostico:", nrow(vals_all), "\n")
  return(vals_all)
}

vals_cedrela      <- extrair_valores_amb(cedrela_pts,      r_stack_raster,
                                         n_background = 10000,
                                         ext = ext_brasil,
                                         especie = "Cedrela")

vals_handroanthus <- extrair_valores_amb(handroanthus_pts, r_stack_raster,
                                         n_background = 10000,
                                         ext = ext_brasil,
                                         especie = "Handroanthus")


#'
#' ## 5.1.2 Teste de multicolinearidade (VIF)
#'
#' Removemos iterativamente a variavel com maior VIF ate que todas tenham
#' VIF <= 5 (limiar adotado por Tesfamariam et al. 2022).
#'
## ----vif----------------------------------------------------------------------
aplicar_vif <- function(vals_amb, limiar_vif = 5, especie) {
  
  cat("\n--- VIF:", especie, "(limiar =", limiar_vif, ") ---\n")
  
  # vifstep: remove iterativamente a variavel com maior VIF ate todas <= limiar
  v <- usdm::vifstep(vals_amb, th = limiar_vif)
  print(v)
  
  vars_mantidas <- v@results$Variables
  cat("Variaveis mantidas apos VIF (", length(vars_mantidas), "):\n",
      paste(vars_mantidas, collapse = ", "), "\n")
  
  return(as.character(vars_mantidas))
}

vars_cedrela_vif      <- aplicar_vif(vals_cedrela,      limiar_vif = 5,
                                     especie = "Cedrela")
vars_handroanthus_vif <- aplicar_vif(vals_handroanthus, limiar_vif = 5,
                                     especie = "Handroanthus")

# Subset do stack mantendo apenas as variaveis aprovadas no VIF
stack_cedrela_vif      <- raster::subset(r_stack_raster, vars_cedrela_vif)
stack_handroanthus_vif <- raster::subset(r_stack_raster, vars_handroanthus_vif)


#'
#' ## 5.1.3 Teste de Jackknife (importancia das variaveis)
#'
#' Apos o VIF, ajustamos um MaxEnt preliminar com `jackknife=true` para
#' identificar variaveis com baixo ganho. Removemos aquelas cuja contribuicao
#' percentual seja inferior ao limiar (ex.: 1%). Ajuste `limiar_contrib`
#' conforme o objetivo do estudo.
#'
## ----jackknife----------------------------------------------------------------
aplicar_jackknife <- function(pontos, stack_vif, especie,
                              limiar_contrib = 1) {
  
  cat("\n--- Jackknife:", especie, "(limiar contribuicao =",
      limiar_contrib, "%) ---\n")
  
  args_maxent <- c(
    "jackknife=true",          # Ativa o teste de Jackknife
    "responsecurves=true",     # Curvas de resposta (auxilia interpretacao)
    "writebackgroundpredictions=false"
  )
  
  xm_prelim <- dismo::maxent(stack_vif, pontos, args = args_maxent)
  
  # Extrair contribuicao percentual de cada variavel
  contrib <- xm_prelim@results
  idx_contrib <- grep("\\.contribution$", rownames(contrib))
  df_contrib <- data.frame(
    variavel    = gsub("\\.contribution$", "", rownames(contrib)[idx_contrib]),
    contribuicao = as.numeric(contrib[idx_contrib, 1])
  )
  df_contrib <- df_contrib[order(-df_contrib$contribuicao), ]
  cat("Contribuicao percentual das variaveis:\n")
  print(df_contrib)
  
  # Variaveis mantidas: contribuicao >= limiar
  vars_mantidas <- df_contrib$variavel[df_contrib$contribuicao >= limiar_contrib]
  cat("\nVariaveis mantidas apos Jackknife (", length(vars_mantidas), "):\n",
      paste(vars_mantidas, collapse = ", "), "\n")
  
  # Plotar grafico de Jackknife do MaxEnt
  plot(xm_prelim, main = paste("Importancia (Jackknife) -", especie))
  
  return(list(
    modelo_prelim = xm_prelim,
    contribuicao  = df_contrib,
    vars_final    = as.character(vars_mantidas)
  ))
}

par(mfrow = c(1, 2))
jk_cedrela      <- aplicar_jackknife(cedrela_pts,      stack_cedrela_vif,
                                     "Cedrela",      limiar_contrib = 1)
jk_handroanthus <- aplicar_jackknife(handroanthus_pts, stack_handroanthus_vif,
                                     "Handroanthus", limiar_contrib = 1)
par(mfrow = c(1, 1))

# Stacks finais (apos VIF + Jackknife)
stack_cedrela_final      <- raster::subset(r_stack_raster, jk_cedrela$vars_final)
stack_handroanthus_final <- raster::subset(r_stack_raster, jk_handroanthus$vars_final)

cat("\n=== Resumo da selecao de variaveis ===\n")
cat("Cedrela     : inicio =", nlyr(r_stack), " -> apos VIF =",
    length(vars_cedrela_vif), " -> apos Jackknife =",
    length(jk_cedrela$vars_final), "\n")
cat("Handroanthus: inicio =", nlyr(r_stack), " -> apos VIF =",
    length(vars_handroanthus_vif), " -> apos Jackknife =",
    length(jk_handroanthus$vars_final), "\n")


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
  set.seed(1350)
  group      <- dismo::kfold(pontos, k)
  pres_train <- pontos[group != 1, ]
  pres_test  <- pontos[group == 1, ]
  
  cat("Treino:", nrow(pres_train), "| Teste:", nrow(pres_test), "\n")
  
  # 6.2 Ajuste MaxEnt (requer maxent.jar)
  cat("Ajustando MaxEnt...\n")
  xm <- dismo::maxent(climate_stack, pres_train)
  
  # 6.3 Geracao de background (pseudo-ausencias)
  set.seed(1350)
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
#' -------------------------------------------------------------------------
#' OPCAO ALTERNATIVA (Tesfamariam et al. 2022):
#' Regressao logistica do MaxEnt executada com 100 replicas e 500 iteracoes
#' para permitir tempo suficiente de convergencia. O mapa de distribuicao e
#' adequabilidade de habitat e produzido pela MEDIA das 100 repeticoes.
#' Para usar esta opcao, descomente o bloco abaixo e troque as chamadas em
#' "7. Ajuste dos modelos" por `ajustar_maxent_replicado`.
#' -------------------------------------------------------------------------
#'
## ----funcao_maxent_replicado--------------------------------------------------
# ajustar_maxent_replicado <- function(pontos, climate_stack, especie,
#                                      ext,
#                                      n_replicas  = 100,
#                                      n_iteracoes = 500,
#                                      pct_teste   = 25) {
#
#   cat("\n========== MaxEnt replicado:", especie, "==========\n")
#   cat("Replicas:", n_replicas, "| Iteracoes:", n_iteracoes,
#       "| % teste:", pct_teste, "\n")
#
#   # Argumentos do MaxEnt seguindo o artigo:
#   #   - replicates=100         -> 100 replicas
#   #   - maximumiterations=500  -> 500 iteracoes por replica
#   #   - replicatetype=subsample -> particao aleatoria treino/teste a cada replica
#   #   - randomtestpoints=25    -> 25% dos pontos para teste, 75% para treino
#   #   - betamultiplier=1       -> regularizacao padrao (reduz overfitting)
#   #   - outputformat=logistic  -> saida em probabilidade logistica
#   #   - jackknife=true         -> mantem o jackknife para diagnostico
#   args_replicado <- c(
#     paste0("replicates=", n_replicas),
#     paste0("maximumiterations=", n_iteracoes),
#     "replicatetype=subsample",
#     paste0("randomtestpoints=", pct_teste),
#     "betamultiplier=1",
#     "outputformat=logistic",
#     "jackknife=true",
#     "responsecurves=true",
#     "writebackgroundpredictions=false"
#   )
#
#   # Ajuste do MaxEnt com replicas
#   xm <- dismo::maxent(climate_stack, pontos, args = args_replicado)
#
#   # Predicao espacial: quando o modelo foi ajustado com replicates > 1,
#   # dismo::predict retorna a MEDIA das replicas como mapa final.
#   cat("Gerando predicao media das", n_replicas, "replicas...\n")
#   p_pred_media <- dismo::predict(climate_stack, xm, ext = ext, progress = "")
#
#   # AUC medio das replicas (extraido de @results)
#   res <- xm@results
#   auc_idx <- grep("Test.AUC", rownames(res))
#   if (length(auc_idx) > 0) {
#     auc_medio <- mean(as.numeric(res[auc_idx, ]), na.rm = TRUE)
#   } else {
#     auc_medio <- NA
#   }
#   cat("AUC medio (teste):", round(auc_medio, 4), "\n")
#
#   return(list(
#     especie    = especie,
#     modelo     = xm,
#     auc_medio  = auc_medio,
#     predicao   = p_pred_media,   # media das 100 replicas
#     n_replicas = n_replicas
#   ))
# }
#
# # ----- Uso da versao replicada (descomente para executar) -----
# # resultado_cedrela <- ajustar_maxent_replicado(
# #   pontos        = cedrela_pts,
# #   climate_stack = stack_cedrela_final,
# #   especie       = "Cedrela",
# #   ext           = ext_brasil,
# #   n_replicas    = 100,
# #   n_iteracoes   = 500,
# #   pct_teste     = 25
# # )
# #
# # resultado_handroanthus <- ajustar_maxent_replicado(
# #   pontos        = handroanthus_pts,
# #   climate_stack = stack_handroanthus_final,
# #   especie       = "Handroanthus",
# #   ext           = ext_brasil,
# #   n_replicas    = 100,
# #   n_iteracoes   = 500,
# #   pct_teste     = 25
# # )


#'
#' # 7. Ajuste dos modelos para cada especie
#'
#' IMPORTANTE: agora usamos os stacks ja reduzidos (apos VIF + Jackknife)
#' `stack_cedrela_final` e `stack_handroanthus_final` em vez do stack completo.
#'
## ----ajuste-cedrela-----------------------------------------------------------
resultado_cedrela <- ajustar_maxent(
  pontos        = cedrela_pts,
  climate_stack = stack_cedrela_final,
  especie       = "Cedrela",
  ext           = ext_brasil,
  n_background  = 1000,
  k             = 5
)

#'
## ----ajuste-handroanthus------------------------------------------------------
resultado_handroanthus <- ajustar_maxent(
  pontos        = handroanthus_pts,
  climate_stack = stack_handroanthus_final,
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
#' ocorrencia sobrepostos e contorno dos estados do Brasil.
#'
## ----mapas-presente-----------------------------------------------------------
# 10.1 Carregar contorno do Brasil (estados) via geobr
brasil_uf <- tryCatch(
  geobr::read_state(year = 2020, showProgress = FALSE),
  error = function(e) {
    warning("Falha ao baixar contorno via geobr: ", conditionMessage(e),
            "\nMapas serao gerados sem contorno dos estados.")
    NULL
  }
)

# Converter para Spatial (compatibilidade com raster::plot + add = TRUE)
if (!is.null(brasil_uf)) {
  brasil_sp <- sf::as_Spatial(sf::st_geometry(brasil_uf))
} else {
  brasil_sp <- NULL
}

# 10.2 Paleta de cores: do azul (baixa adequabilidade) ao verde (alta),
# similar a Tesfamariam et al. (2022, Fig. 2a)
paleta_adeq <- colorRampPalette(c(
  "#08306b", "#2171b5", "#6baed6", "#c6dbef",
  "#ffffcc", "#a1d99b", "#41ab5d", "#005a32"
))(100)

# 10.3 Funcao para plotar o mapa de adequabilidade de uma especie
plot_mapa_adeq <- function(predicao, pontos, especie,
                           contorno = NULL, paleta = paleta_adeq) {
  
  raster::plot(predicao,
               col   = paleta,
               main  = paste("Adequabilidade de habitat -", especie),
               xlab  = "Longitude",
               ylab  = "Latitude",
               legend.args = list(text = "Adequabilidade",
                                  side = 4, line = 2.5, cex = 0.8))
  
  # Contorno dos estados (se disponivel)
  if (!is.null(contorno)) {
    sp::plot(contorno, add = TRUE, border = "gray30", lwd = 0.5)
  }
  
  # Pontos de ocorrencia
  points(pontos$lon, pontos$lat,
         pch = 21, bg = "red", col = "black", cex = 0.6)
  
  # Legenda dos pontos
  legend("bottomright",
         legend = "Ocorrencias",
         pch = 21, pt.bg = "red", col = "black",
         bty = "n", cex = 0.8)
}

# 10.4 Mapa Cedrela
plot_mapa_adeq(predicao = resultado_cedrela$predicao,
               pontos   = cedrela_pts,
               especie  = "Cedrela",
               contorno = brasil_sp)

# 10.5 Mapa Handroanthus
plot_mapa_adeq(predicao = resultado_handroanthus$predicao,
               pontos   = handroanthus_pts,
               especie  = "Handroanthus",
               contorno = brasil_sp)

# 10.6 (Opcional) Mapa binario: adequado vs nao-adequado
#
# Seguindo o artigo, classificamos como adequado todo pixel com adequabilidade
# >= MENOR valor predito nos pontos de presenca (criterio "minimum training
# presence" / "minimum predicted area"). Pixels abaixo desse limiar viram
# nao-adequados.
#
classificar_binario <- function(predicao, pontos, especie) {
  
  # Adequabilidade nos pontos de presenca
  vals_pres <- raster::extract(predicao, pontos, method="bilinear")
  vals_pres <- vals_pres[!is.na(vals_pres)]
  
  limiar <- min(vals_pres)
  cat("Limiar (min training presence) para", especie, ":",
      round(limiar, 4), "\n")
  
  bin <- raster::reclassify(predicao,
                            rcl = matrix(c(-Inf, limiar, 0,
                                           limiar,  Inf, 1),
                                         ncol = 3, byrow = TRUE))
  return(list(binario = bin, limiar = limiar))
}

bin_cedrela      <- classificar_binario(resultado_cedrela$predicao,
                                        cedrela_pts,      "Cedrela")
bin_handroanthus <- classificar_binario(resultado_handroanthus$predicao,
                                        handroanthus_pts, "Handroanthus")

# Paleta binaria: cinza claro (nao-adequado) e verde escuro (adequado)
paleta_bin <- c("#d9d9d9", "#005a32")

par(mfrow = c(1, 2))
raster::plot(bin_cedrela$binario,
             col = paleta_bin,
             main = "Cedrela - adequado vs nao-adequado",
             legend = FALSE, xlab = "Longitude", ylab = "Latitude")
if (!is.null(brasil_sp)) sp::plot(brasil_sp, add = TRUE,
                                  border = "gray30", lwd = 0.5)
points(cedrela_pts$lon, cedrela_pts$lat,
       pch = 21, bg = "red", col = "black", cex = 0.5)
legend("bottomright",
       legend = c("Nao-adequado", "Adequado", "Ocorrencias"),
       fill   = c(paleta_bin, NA),
       pch    = c(NA, NA, 21),
       pt.bg  = c(NA, NA, "red"),
       border = c("black", "black", NA),
       bty = "n", cex = 0.7)

raster::plot(bin_handroanthus$binario,
             col = paleta_bin,
             main = "Handroanthus - adequado vs nao-adequado",
             legend = FALSE, xlab = "Longitude", ylab = "Latitude")
if (!is.null(brasil_sp)) sp::plot(brasil_sp, add = TRUE,
                                  border = "gray30", lwd = 0.5)
points(handroanthus_pts$lon, handroanthus_pts$lat,
       pch = 21, bg = "red", col = "black", cex = 0.5)
legend("bottomright",
       legend = c("Nao-adequado", "Adequado", "Ocorrencias"),
       fill   = c(paleta_bin, NA),
       pch    = c(NA, NA, 21),
       pt.bg  = c(NA, NA, "red"),
       border = c("black", "black", NA),
       bty = "n", cex = 0.7)
par(mfrow = c(1, 1))

# 10.7 Percentual da area classificada como adequada (similar ao "> 48%"
# reportado no artigo)
calc_pct_adequado <- function(bin_raster, especie) {
  freqs <- raster::freq(bin_raster, useNA = "no")
  total <- sum(freqs[, "count"])
  adeq  <- freqs[freqs[, "value"] == 1, "count"]
  pct   <- 100 * adeq / total
  cat(especie, ": area adequada =", round(pct, 2), "%\n")
  return(pct)
}

pct_cedrela      <- calc_pct_adequado(bin_cedrela$binario,      "Cedrela")
pct_handroanthus <- calc_pct_adequado(bin_handroanthus$binario, "Handroanthus")


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

# Exportar rasters de predicao (continuos)
raster::writeRaster(resultado_cedrela$predicao,
                    filename = file.path(out_dir, "maxent_cedrela_presente.tif"),
                    overwrite = TRUE)

raster::writeRaster(resultado_handroanthus$predicao,
                    filename = file.path(out_dir, "maxent_handroanthus_presente.tif"),
                    overwrite = TRUE)

# Exportar rasters binarios (adequado vs nao-adequado)
raster::writeRaster(bin_cedrela$binario,
                    filename = file.path(out_dir, "maxent_cedrela_binario.tif"),
                    overwrite = TRUE)

raster::writeRaster(bin_handroanthus$binario,
                    filename = file.path(out_dir, "maxent_handroanthus_binario.tif"),
                    overwrite = TRUE)

# Salvar modelos R (para reabrir sem refazer o ajuste)
saveRDS(resultado_cedrela,      file.path(out_dir, "modelo_cedrela.rds"))
saveRDS(resultado_handroanthus, file.path(out_dir, "modelo_handroanthus.rds"))

# Tabela-resumo dos AUC + variaveis selecionadas + % de area adequada
resumo <- data.frame(
  Especie  = c("Cedrela", "Handroanthus"),
  AUC      = c(resultado_cedrela$auc, resultado_handroanthus$auc),
  N_treino = c(nrow(resultado_cedrela$pres_train),
               nrow(resultado_handroanthus$pres_train)),
  N_teste  = c(nrow(resultado_cedrela$pres_test),
               nrow(resultado_handroanthus$pres_test)),
  Vars_apos_VIF      = c(length(vars_cedrela_vif),
                         length(vars_handroanthus_vif)),
  Vars_apos_Jacknife = c(length(jk_cedrela$vars_final),
                         length(jk_handroanthus$vars_final)),
  Vars_finais        = c(paste(jk_cedrela$vars_final,      collapse = "; "),
                         paste(jk_handroanthus$vars_final, collapse = "; ")),
  Limiar_binario     = c(round(bin_cedrela$limiar,      4),
                         round(bin_handroanthus$limiar, 4)),
  Pct_area_adequada  = c(round(pct_cedrela,      2),
                         round(pct_handroanthus, 2))
)
print(resumo)

writexl::write_xlsx(resumo, file.path(out_dir, "resumo_maxent.xlsx"))