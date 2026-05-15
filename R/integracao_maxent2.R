# =========================================================================
# Modelos de nicho ecologico para Cedrela e Handroanthus no Brasil (MaxEnt)
#
# Pipeline ajustado para reproduzir a metodologia de:
#   Elera-Gonzales, D.G. et al. (2026). Potential geographic distribution of
#   Handroanthus serratifolius in tropical America. Flora 336, 152926.
#   https://doi.org/10.1016/j.flora.2026.152926
#
# Mantemos: mesma area de estudo (Brasil), mesmos rasters de `rasters/`,
#           mesmos dados (cedrela_br_var_amb.csv, handroanthus_var_amb.csv).
#
# Adaptamos do artigo:
#   - Thinning espacial de 5 km nas ocorrencias (spThin)  -- Sec. 2.1
#   - Selecao de variaveis por Pearson |r| < 0.70 + VIF < 10  -- Sec. 2.3
#   - MaxEnt: 100 replicas, 5000 iteracoes, Bootstrap, 75% treino / 25% teste,
#     saida Cloglog, regra "Maximum training sensitivity plus specificity"
#     (maxSSS)  -- Sec. 2.3
#   - Avaliacao por AUC treino/teste (com SD das 100 replicas) e TSS no maxSSS
#     -- Sec. 2.4
#   - Reclassificacao em 4 classes (limiares 0.265, 0.53, 0.765) e mapa
#     binario no maxSSS  -- Sec. 2.4 / 3.3
#   - Calculo de area em km^2 (raster::area corrige pela latitude)
# =========================================================================

#'
## ----setup--------------------------------------------------------------------
#| include: false

# Ajuste global para graficos
library(knitr)
opts_knit$set(global.par = TRUE)
par(mar = c(5, 5, 1, 1))

# ---------------------------------------------------------------------------
# 1. Configuracao inicial
# ---------------------------------------------------------------------------

# Limpar ambiente (descomente se desejar)
# rm(list = ls()); gc(reset = TRUE); graphics.off()

# Heap do Java -- MaxEnt replicado consome bastante memoria
options(java.parameters = "-Xmx8g")

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
if(!require(spThin))      install.packages("spThin",      dep = TRUE, quiet = TRUE)
# 'usdm' fornece vifcor()/vifstep() para diagnostico de colinearidade
# 'spThin' aplica filtro espacial baseado em distancia minima entre pontos

# Verificacao do binario maxent.jar
maxent_jar <- file.path(system.file("java", package = "dismo"), "maxent.jar")
if (!file.exists(maxent_jar)) {
  warning("maxent.jar nao encontrado em: ", maxent_jar,
          "\nBaixe em https://biodiversityinformatics.amnh.org/open_source/maxent/")
}


#'
#' # 2. Reprodutibilidade
#'
## ----seed---------------------------------------------------------------------
set.seed(1350)


#'
#' # 3. Carregando os conjuntos de dados
#'
## ----importacao_dados---------------------------------------------------------
cedrela      <- read.csv(here::here("data", "cedrela_br_var_amb.csv"))
handroanthus <- read.csv(here::here("data", "handroanthus_var_amb.csv"))
head(cedrela); head(handroanthus)


#'
#' # 4. Carregamento do stack de rasters
#'
## ----stack--------------------------------------------------------------------
# raster_files <- list.files(here::here("rasters"),
#                            pattern = "\\.tif$", full.names = TRUE)

raster_files <- list.files(here::here("rasters"),   pattern = "32\\.vapor|33\\.solar\\.tif$",   full.names = TRUE)

# Excluir stacks consolidados e recortes regionais (mesmas exclusoes do original)
excluir_pattern <- "Brazil_env_stack|Amazon_stack|Caatinga_elev_corridor|elev_corridor_minus1|48.density"
raster_files <- raster_files[!grepl(excluir_pattern, raster_files)]

# Empilhar como SpatRaster (inspecao) e como RasterStack (dismo/MaxEnt)
r_stack <- terra::rast(raster_files)
names(r_stack) <- tools::file_path_sans_ext(basename(raster_files))

r_stack_raster <- raster::stack(raster_files)
names(r_stack_raster) <- tools::file_path_sans_ext(basename(raster_files))

cat("Camadas no stack:", nlyr(r_stack), "\n")
cat("Extensao espacial:\n"); print(terra::ext(r_stack))

# Extensao do estudo (Brasil)
ext_brasil <- raster::extent(r_stack_raster)


#'
#' # 5. Preparacao dos dados de ocorrencia
#'
#' Para cada especie:
#'   (i)  extrair coordenadas (lon, lat)
#'   (ii) remover NAs, duplicatas e coordenadas implausiveis
#'   (iii) aplicar filtro espacial de 5 km (spThin) - SEC. 2.1 DO ARTIGO
#'
## ----preparacao---------------------------------------------------------------
preparar_ocorrencias <- function(dados, especie) {
  
  set.seed(1350)
  
  cat("\n--- Preparando:", especie, "---\n")
  
  # Padronizar nomes das colunas de coordenadas
  if (all(c("x", "y") %in% colnames(dados))) {
    pontos <- data.frame(lon = dados$x, lat = dados$y)
  } else if (all(c("longitude", "latitude") %in% colnames(dados))) {
    pontos <- data.frame(lon = dados$longitude, lat = dados$latitude)
  } else {
    stop("Colunas de coordenadas nao encontradas para ", especie)
  }
  
  # Remover NAs, duplicatas e coordenadas absurdas
  pontos <- pontos[complete.cases(pontos), ]
  pontos <- pontos[!duplicated(pontos), ]
  pontos <- pontos[pontos$lon >= -180 & pontos$lon <= 180 &
                     pontos$lat >=  -90 & pontos$lat <=  90, ]
  
  cat("Pontos validos (antes do thinning):", nrow(pontos), "\n")
  return(pontos)
}


#'
#' ## 5.1 Filtro espacial de 5 km (Sec. 2.1 do artigo)
#'
#' O artigo aplica esse filtro com a ferramenta "Near" do ArcGIS Pro. No R,
#' usamos spThin (Aiello-Lammens et al., 2015) que mantem apenas pontos
#' separados por uma distancia minima especificada.
#'
## ----thinning-----------------------------------------------------------------
aplicar_thinning <- function(pontos, especie, dist_km = 5,
                             n_reps = 50, seed = 1350) {
  
  cat("--- Thinning espacial de", dist_km, "km:", especie, "---\n")
  pts <- pontos
  pts$especie <- especie
  
  set.seed(seed)
  thin_out <- spThin::thin(
    loc.data           = pts,
    lat.col            = "lat",
    long.col           = "lon",
    spec.col           = "especie",
    thin.par           = dist_km,
    reps               = n_reps,
    locs.thinned.list.return = TRUE,
    write.files        = FALSE,
    write.log.file     = FALSE,
    verbose            = FALSE
  )
  
  # spThin retorna varios subconjuntos -- mantemos o maior (max pontos retidos)
  idx_maior <- which.max(sapply(thin_out, nrow))
  out <- thin_out[[idx_maior]]
  names(out) <- c("lon", "lat")
  
  cat("Pontos retidos apos thinning:", nrow(out), "\n")
  return(out)
}

cedrela_raw      <- preparar_ocorrencias(cedrela,      "Cedrela")
handroanthus_raw <- preparar_ocorrencias(handroanthus, "Handroanthus")

cedrela_pts      <- aplicar_thinning(cedrela_raw,      "Cedrela",      5)
handroanthus_pts <- aplicar_thinning(handroanthus_raw, "Handroanthus", 5)

# Visualizar distribuicao geografica
par(mfrow = c(1, 2))
plot(cedrela_pts$lon, cedrela_pts$lat,
     col = "darkgreen", pch = 19, cex = 0.5,
     xlab = "Longitude", ylab = "Latitude",
     main = paste0("Cedrela (n=", nrow(cedrela_pts), ")"))
maps::map(add = TRUE)

plot(handroanthus_pts$lon, handroanthus_pts$lat,
     col = "darkorange", pch = 19, cex = 0.5,
     xlab = "Longitude", ylab = "Latitude",
     main = paste0("Handroanthus (n=", nrow(handroanthus_pts), ")"))
maps::map(add = TRUE)
par(mfrow = c(1, 1))


#'
#' # 6. Selecao de variaveis ambientais: Pearson |r| < 0.70 + VIF < 10
#'
#' SUBSTITUI o passo VIF<5 + Jackknife do script original.
#' Segue o procedimento da Sec. 2.3 do artigo (Paul et al., 2021; Naimi, 2023):
#'   (a) usdm::vifcor() -- remove iterativamente variaveis com |r| > 0.70
#'   (b) usdm::vifstep() sobre o conjunto resultante -- garante VIF <= 10
#'
## ----extrair-valores----------------------------------------------------------
extrair_valores_amb <- function(pontos, climate_stack, n_background = 10000,
                                ext, especie, seed = 1350) {
  
  cat("\n--- Extraindo valores ambientais:", especie, "---\n")
  
  set.seed(seed)
  backg <- dismo::randomPoints(climate_stack, n = n_background,
                               ext = ext, extf = 1)
  colnames(backg) <- c("lon", "lat")
  
  vals_pres <- raster::extract(climate_stack, pontos[, c("lon", "lat")],
                               method = "bilinear")
  vals_back <- raster::extract(climate_stack, backg, method = "bilinear")
  
  vals_all <- as.data.frame(rbind(vals_pres, vals_back))
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


## ----selecao-vars-------------------------------------------------------------
selecionar_variaveis <- function(vals_amb, especie,
                                 cor_th = 0.70, vif_th = 5) {
  
  cat("\n--- Selecao de variaveis:", especie, " (|r|<", cor_th,
      "e VIF<", vif_th, ") ---\n")
  
  # (a) Pearson |r| < 0.70  -- vifcor
  sel_cor <- usdm::vifcor(vals_amb, th = cor_th, method="spearman")
  vars_cor <- as.character(sel_cor@results$Variables)
  cat("Apos correlacao Pearson < ", cor_th, ":", length(vars_cor), "variaveis\n")
  
  # (b) VIF < 10  -- vifstep
  sel_vif <- usdm::vifstep(vals_amb[, vars_cor], th = vif_th)
  vars_final <- as.character(sel_vif@results$Variables)
  cat("Apos VIF <", vif_th, ":", length(vars_final), "variaveis\n")
  cat("Variaveis finais:\n", paste(vars_final, collapse = ", "), "\n")
  
  return(list(
    vars_cor      = vars_cor,
    vars_final    = vars_final,
    info_vifcor   = sel_cor,
    info_vifstep  = sel_vif
  ))
}

sel_cedrela      <- selecionar_variaveis(vals_cedrela,      "Cedrela")
sel_handroanthus <- selecionar_variaveis(vals_handroanthus, "Handroanthus")

# Stacks finais
stack_cedrela_final      <- raster::subset(r_stack_raster, sel_cedrela$vars_final)
stack_handroanthus_final <- raster::subset(r_stack_raster, sel_handroanthus$vars_final)

cat("\n=== Resumo da selecao de variaveis ===\n")
cat("Cedrela     : inicio =", nlyr(r_stack),
    " -> apos |r|<0.70 =", length(sel_cedrela$vars_cor),
    " -> apos VIF<10 =",   length(sel_cedrela$vars_final), "\n")
cat("Handroanthus: inicio =", nlyr(r_stack),
    " -> apos |r|<0.70 =", length(sel_handroanthus$vars_cor),
    " -> apos VIF<10 =",   length(sel_handroanthus$vars_final), "\n")


#'
#' # 7. Funcao para ajuste do MaxEnt REPLICADO (configuracao do artigo)
#'
#' Configuracao identica a do artigo (Sec. 2.3):
#'   - replicates       = 100        100 repeticoes
#'   - maximumiterations = 5000      5000 iteracoes por replica
#'   - replicatetype    = bootstrap  amostragem com reposicao
#'   - randomtestpoints = 25         25% para teste, 75% para treino
#'   - outputformat     = cloglog    Cloglog (recomendado >= MaxEnt 3.4.0)
#'   - threshold rule   = Maximum training Sensitivity Plus Specificity (maxSSS)
#'   - jackknife        = true       teste de Jackknife (importancia de variaveis)
#'   - responsecurves   = true       curvas de resposta
#'   - betamultiplier   = 1          regularizacao padrao
#'
#' Quando replicates > 1, dismo::predict retorna automaticamente o MAPA MEDIO
#' das replicas, exatamente como descrito no artigo.
#'
## ----funcao-maxent-replicado--------------------------------------------------
ajustar_maxent <- function(pontos, climate_stack, especie, ext,
                           dir_saida_modelo,
                           n_replicas  = 10,
                           n_iteracoes = 50,
                           pct_teste   = 25) {
  
  cat("\n========== MaxEnt:", especie, "==========\n")
  cat("Replicas:", n_replicas, "| Iteracoes:", n_iteracoes,
      "| % teste:", pct_teste, "\n")
  
  if (!dir.exists(dir_saida_modelo)) {
    dir.create(dir_saida_modelo, recursive = TRUE)
  }
  
  args_replicado <- c(
    paste0("replicates=", n_replicas),
    paste0("maximumiterations=", n_iteracoes),
    "replicatetype=bootstrap",
    paste0("randomtestpoints=", pct_teste),
    "outputformat=cloglog",
    "betamultiplier=1",
    "threshold=true",
    "jackknife=true",
    "responsecurves=true",
    "writeplotdata=true",
    "writebackgroundpredictions=false",
    "applythresholdrule=Maximum training sensitivity plus specificity"
  )
  
  pts_xy <- pontos[, c("lon", "lat")]
  
  cat("Ajustando MaxEnt (pode demorar varios minutos)...\n")
  xm <- dismo::maxent(climate_stack, pts_xy,
                      args = args_replicado,
                      path = dir_saida_modelo,
                      silent = TRUE)
  
  # ---- Metricas medias das replicas ----
  res <- xm@results
  
  pegar_metrica <- function(nome) {
    idx <- grep(nome, rownames(res), fixed = FALSE)
    if (length(idx) == 0) return(c(media = NA, sd = NA))
    vals <- as.numeric(res[idx, ])
    vals <- vals[!is.na(vals)]
    n_use <- min(n_replicas, length(vals))
    c(media = mean(vals[1:n_use], na.rm = TRUE),
      sd    = sd(vals[1:n_use],   na.rm = TRUE))
  }
  
  auc_treino <- pegar_metrica("Training.AUC")
  auc_teste  <- pegar_metrica("Test.AUC")
  auc_geral  <- mean(c(auc_treino["media"], auc_teste["media"]), na.rm = TRUE)
  
  cat("\n--- Performance media (", n_replicas, " replicas) ---\n")
  cat("AUC treino: ", round(auc_treino["media"], 3),
      " (SD = ", round(auc_treino["sd"], 3), ")\n", sep = "")
  cat("AUC teste:  ", round(auc_teste["media"],  3),
      " (SD = ", round(auc_teste["sd"],  3), ")\n", sep = "")
  cat("AUC geral:  ", round(auc_geral, 3), "\n", sep = "")
  
  # ---- Predicao espacial (mapa medio das replicas) ----
  cat("\nGerando mapa medio de adequabilidade...\n")
  p_pred <- dismo::predict(climate_stack, xm, ext = ext, progress = "text")
  
  return(list(
    especie     = especie,
    modelo      = xm,
    predicao    = p_pred,
    auc_treino  = auc_treino,
    auc_teste   = auc_teste,
    auc_geral   = auc_geral,
    n_replicas  = n_replicas,
    resultados  = res
  ))
}


#'
#' # 8. Ajuste dos modelos para cada especie
#'
## ----ajuste-cedrela-----------------------------------------------------------
out_dir <- here::here("maxent_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

resultado_cedrela <- ajustar_maxent(
  pontos          = cedrela_pts,
  climate_stack   = stack_cedrela_final,
  especie         = "Cedrela",
  ext             = ext_brasil,
  dir_saida_modelo = file.path(out_dir, "maxent_cedrela"),
  n_replicas      = 1,
  n_iteracoes     = 5,
  pct_teste       = 25
)

## ----ajuste-handroanthus------------------------------------------------------
resultado_handroanthus <- ajustar_maxent(
  pontos          = handroanthus_pts,
  climate_stack   = stack_handroanthus_final,
  especie         = "Handroanthus",
  ext             = ext_brasil,
  dir_saida_modelo = file.path(out_dir, "maxent_handroanthus"),
  n_replicas      = 1,
  n_iteracoes     = 5,
  pct_teste       = 25
)


#'
#' # 9. Avaliacao adicional: TSS no limiar maxSSS (Liu et al., 2016)
#'
#' TSS = sensibilidade + especificidade - 1  no limiar Cloglog correspondente a
#' "Maximum training Sensitivity Plus Specificity" (maxSSS). Esse e o mesmo
#' criterio adotado pelo artigo (Sec. 2.4), que reportou TSS = 0.436 ao limiar
#' P = 0.53. Para nossas especies e dados, o limiar e o TSS serao especificos.
#'
## ----tss----------------------------------------------------------------------
calcular_tss <- function(modelo, predicao, pontos, climate_stack, ext,
                         n_bg = 10000, n_replicas = 100, seed = 1350) {
  
  res <- modelo@results
  
  # Limiar maxSSS Cloglog -- media das replicas
  idx_thr <- grep("Maximum.training.sensitivity.plus.specificity.Cloglog.threshold",
                  rownames(res), fixed = FALSE)
  if (length(idx_thr) == 0) {
    idx_thr <- grep("Maximum.training.sensitivity.plus.specificity",
                    rownames(res), fixed = FALSE)
  }
  thr_vals <- as.numeric(res[idx_thr, ])
  thr_vals <- thr_vals[!is.na(thr_vals)]
  thr_maxSSS <- mean(thr_vals[1:min(n_replicas, length(thr_vals))])
  
  cat("Limiar maxSSS (media replicas):", round(thr_maxSSS, 3), "\n")
  
  # Background para avaliacao
  set.seed(seed)
  bg <- dismo::randomPoints(climate_stack, n = n_bg, ext = ext)
  colnames(bg) <- c("lon", "lat")
  
  s_pres <- raster::extract(predicao, pontos[, c("lon", "lat")])
  s_bg   <- raster::extract(predicao, bg)
  s_pres <- s_pres[!is.na(s_pres)]
  s_bg   <- s_bg[!is.na(s_bg)]
  
  # Matriz de confusao no limiar maxSSS
  TP <- sum(s_pres >= thr_maxSSS)
  FN <- sum(s_pres <  thr_maxSSS)
  FP <- sum(s_bg   >= thr_maxSSS)
  TN <- sum(s_bg   <  thr_maxSSS)
  
  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  TSS  <- sens + spec - 1
  om   <- FN / (TP + FN)
  
  cat("Sensibilidade:", round(sens, 3),
      "| Especificidade:", round(spec, 3),
      "| TSS:", round(TSS, 3),
      "| Omissao:", round(om, 3), "\n")
  
  return(list(
    limiar        = thr_maxSSS,
    sensibilidade = sens,
    especificidade = spec,
    TSS           = TSS,
    omissao       = om
  ))
}

cat("\n--- TSS Cedrela ---\n")
tss_cedrela <- calcular_tss(resultado_cedrela$modelo,
                            resultado_cedrela$predicao,
                            cedrela_pts,
                            stack_cedrela_final,
                            ext_brasil)

cat("\n--- TSS Handroanthus ---\n")
tss_handroanthus <- calcular_tss(resultado_handroanthus$modelo,
                                 resultado_handroanthus$predicao,
                                 handroanthus_pts,
                                 stack_handroanthus_final,
                                 ext_brasil)


#'
#' # 10. Importancia das variaveis e curvas de resposta
#'
## ----importancia--------------------------------------------------------------
par(mfrow = c(1, 2))
plot(resultado_cedrela$modelo,
     main = "Importancia - Cedrela")
plot(resultado_handroanthus$modelo,
     main = "Importancia - Handroanthus")
par(mfrow = c(1, 1))

# As curvas de resposta e plots de Jackknife sao gerados pelo MaxEnt em:
#   maxent_output/maxent_cedrela/plots/
#   maxent_output/maxent_handroanthus/plots/


#'
#' # 11. Reclassificacao em 4 classes + mapa binario no maxSSS
#'
#' Esquema do artigo (Sec. 2.4):
#'   Alta      (P >= 0.765)
#'   Media     (0.53  <= P < 0.765)
#'   Baixa     (0.265 <= P < 0.53)
#'   Inadeq.   (P < 0.265)
#'
#' Mantemos os limiares fixos do artigo por padrao. Caso prefira limiares
#' especificos por especie, use os valores de tss_*$limiar e derive
#' lim_inf = limiar/2, lim_alto = (limiar+1)/2.
#'
## ----reclassificacao----------------------------------------------------------
reclassificar_4_classes <- function(pred,
                                    lim_inf  = 0.265,
                                    lim_med  = 0.53,
                                    lim_alto = 0.765) {
  rcl <- matrix(c(
    -Inf,    lim_inf,  0,   # inadequado
    lim_inf, lim_med,  1,   # baixa
    lim_med, lim_alto, 2,   # media
    lim_alto, Inf,     3    # alta
  ), ncol = 3, byrow = TRUE)
  raster::reclassify(pred, rcl = rcl)
}

mapa4_cedrela      <- reclassificar_4_classes(resultado_cedrela$predicao)
mapa4_handroanthus <- reclassificar_4_classes(resultado_handroanthus$predicao)

# Mapa binario no limiar maxSSS especifico de cada especie
mapa_bin_cedrela <- raster::reclassify(
  resultado_cedrela$predicao,
  rcl = matrix(c(-Inf, tss_cedrela$limiar, 0,
                 tss_cedrela$limiar, Inf,  1),
               ncol = 3, byrow = TRUE)
)
mapa_bin_handroanthus <- raster::reclassify(
  resultado_handroanthus$predicao,
  rcl = matrix(c(-Inf, tss_handroanthus$limiar, 0,
                 tss_handroanthus$limiar, Inf,  1),
               ncol = 3, byrow = TRUE)
)


#'
#' # 12. Calculo de area por classe (em km^2, corrigida pela latitude)
#'
## ----areas--------------------------------------------------------------------
calcular_area_classes <- function(mapa, nomes_classes) {
  area_pix <- raster::area(mapa, na.rm = TRUE)
  vals     <- raster::values(mapa)
  area_v   <- raster::values(area_pix)
  ok <- !is.na(vals) & !is.na(area_v)
  vals <- vals[ok]; area_v <- area_v[ok]
  
  ag <- tapply(area_v, vals, sum)
  total <- sum(ag)
  
  data.frame(
    classe   = nomes_classes[as.integer(names(ag)) + 1],
    valor    = as.integer(names(ag)),
    area_km2 = round(as.numeric(ag), 0),
    pct      = round(100 * as.numeric(ag) / total, 2),
    stringsAsFactors = FALSE
  )
}

nomes_4cl <- c("Inadequado", "Baixa", "Media", "Alta")
nomes_bin <- c("Inadequado", "Adequado")

areas4_cedrela      <- calcular_area_classes(mapa4_cedrela,      nomes_4cl)
areas4_handroanthus <- calcular_area_classes(mapa4_handroanthus, nomes_4cl)
areasbin_cedrela    <- calcular_area_classes(mapa_bin_cedrela,    nomes_bin)
areasbin_handroanthus <- calcular_area_classes(mapa_bin_handroanthus, nomes_bin)

cat("\n--- Areas Cedrela (4 classes) ---\n");      print(areas4_cedrela)
cat("\n--- Areas Cedrela (binario) ---\n");        print(areasbin_cedrela)
cat("\n--- Areas Handroanthus (4 classes) ---\n"); print(areas4_handroanthus)
cat("\n--- Areas Handroanthus (binario) ---\n");   print(areasbin_handroanthus)


#'
#' # 13. Mapas (contorno do Brasil via geobr, como no original)
#'
## ----contorno-----------------------------------------------------------------
brasil_uf <- tryCatch(
  geobr::read_state(year = 2020, showProgress = FALSE),
  error = function(e) {
    warning("Falha ao baixar contorno via geobr: ", conditionMessage(e),
            "\nMapas serao gerados sem contorno dos estados.")
    NULL
  }
)
brasil_sp <- if (!is.null(brasil_uf))
  sf::as_Spatial(sf::st_geometry(brasil_uf)) else NULL


## ----paletas------------------------------------------------------------------
# Paleta continua (similar a paleta do artigo: azul-verde)
paleta_cont <- colorRampPalette(c(
  "#08306b", "#2171b5", "#6baed6", "#c6dbef",
  "#ffffcc", "#a1d99b", "#41ab5d", "#005a32"
))(100)

# Paleta 4 classes
cores_4cl <- c("#f0f0f0",  # Inadequado
               "#fee08b",  # Baixa
               "#a6d96a",  # Media
               "#1a9850")  # Alta

# Paleta binaria
cores_bin <- c("#d9d9d9", "#c2185b")


## ----plot-continuo------------------------------------------------------------
plot_mapa_continuo <- function(predicao, pontos, especie,
                               contorno = NULL, paleta = paleta_cont) {
  raster::plot(predicao,
               col   = paleta,
               main  = paste("Adequabilidade - ", especie),
               xlab  = "Longitude", ylab = "Latitude",
               zlim  = c(0, 1),
               legend.args = list(text = "Adequabilidade",
                                  side = 4, line = 2.5, cex = 0.8))
  if (!is.null(contorno)) sp::plot(contorno, add = TRUE,
                                   border = "gray30", lwd = 0.5)
  points(pontos$lon, pontos$lat,
         pch = 21, bg = "red", col = "black", cex = 0.5)
  legend("bottomright", legend = "Ocorrencias",
         pch = 21, pt.bg = "red", col = "black",
         bty = "n", cex = 0.8)
}

# Mapas continuos
plot_mapa_continuo(resultado_cedrela$predicao,
                   cedrela_pts, "Cedrela", brasil_sp)
plot_mapa_continuo(resultado_handroanthus$predicao,
                   handroanthus_pts, "Handroanthus", brasil_sp)


## ----plot-4classes------------------------------------------------------------
plot_mapa_4classes <- function(mapa4, especie, contorno = NULL) {
  raster::plot(mapa4,
               col    = cores_4cl,
               breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5),
               main   = paste("Classes de adequabilidade -", especie),
               xlab   = "Longitude", ylab = "Latitude",
               legend = FALSE)
  if (!is.null(contorno)) sp::plot(contorno, add = TRUE,
                                   border = "gray30", lwd = 0.5)
  legend("bottomright",
         legend = c("Inadequado (<0.265)",
                    "Baixa (0.265-0.53)",
                    "Media (0.53-0.765)",
                    "Alta (>=0.765)"),
         fill   = cores_4cl, bty = "n", cex = 0.75)
}

par(mfrow = c(1, 2))
plot_mapa_4classes(mapa4_cedrela,      "Cedrela",      brasil_sp)
plot_mapa_4classes(mapa4_handroanthus, "Handroanthus", brasil_sp)
par(mfrow = c(1, 1))


## ----plot-binario-------------------------------------------------------------
plot_mapa_binario <- function(mapa_bin, pontos, especie,
                              limiar, contorno = NULL) {
  raster::plot(mapa_bin,
               col    = cores_bin,
               breaks = c(-0.5, 0.5, 1.5),
               main   = sprintf("%s - binario (maxSSS = %.3f)",
                                especie, limiar),
               xlab   = "Longitude", ylab = "Latitude",
               legend = FALSE)
  if (!is.null(contorno)) sp::plot(contorno, add = TRUE,
                                   border = "gray30", lwd = 0.5)
  points(pontos$lon, pontos$lat,
         pch = 21, bg = "yellow", col = "black", cex = 0.4)
  legend("bottomright",
         legend = c("Inadequado", "Adequado", "Ocorrencias"),
         fill   = c(cores_bin, NA),
         pch    = c(NA, NA, 21),
         pt.bg  = c(NA, NA, "yellow"),
         border = c("black", "black", NA),
         bty = "n", cex = 0.7)
}

par(mfrow = c(1, 2))
plot_mapa_binario(mapa_bin_cedrela, cedrela_pts, "Cedrela",
                  tss_cedrela$limiar, brasil_sp)
plot_mapa_binario(mapa_bin_handroanthus, handroanthus_pts, "Handroanthus",
                  tss_handroanthus$limiar, brasil_sp)
par(mfrow = c(1, 1))


#'
#' # 14. Exportacao de rasters, modelos e tabelas-resumo
#'
## ----exportacao---------------------------------------------------------------
# Rasters continuos
raster::writeRaster(resultado_cedrela$predicao,
                    file.path(out_dir, "cedrela_adequabilidade_continua.tif"),
                    overwrite = TRUE)
raster::writeRaster(resultado_handroanthus$predicao,
                    file.path(out_dir, "handroanthus_adequabilidade_continua.tif"),
                    overwrite = TRUE)

# Rasters em 4 classes
raster::writeRaster(mapa4_cedrela,
                    file.path(out_dir, "cedrela_4_classes.tif"),
                    overwrite = TRUE, datatype = "INT1U")
raster::writeRaster(mapa4_handroanthus,
                    file.path(out_dir, "handroanthus_4_classes.tif"),
                    overwrite = TRUE, datatype = "INT1U")

# Rasters binarios
raster::writeRaster(mapa_bin_cedrela,
                    file.path(out_dir, "cedrela_binario_maxSSS.tif"),
                    overwrite = TRUE, datatype = "INT1U")
raster::writeRaster(mapa_bin_handroanthus,
                    file.path(out_dir, "handroanthus_binario_maxSSS.tif"),
                    overwrite = TRUE, datatype = "INT1U")

# Modelos R
saveRDS(resultado_cedrela,      file.path(out_dir, "modelo_cedrela.rds"))
saveRDS(resultado_handroanthus, file.path(out_dir, "modelo_handroanthus.rds"))
saveRDS(tss_cedrela,            file.path(out_dir, "tss_cedrela.rds"))
saveRDS(tss_handroanthus,       file.path(out_dir, "tss_handroanthus.rds"))

# Tabela-resumo (estilo Tabela 1/2 do artigo)
resumo <- data.frame(
  Especie  = c("Cedrela", "Handroanthus"),
  N_pontos_pos_thinning_5km = c(nrow(cedrela_pts), nrow(handroanthus_pts)),
  N_vars_pos_cor            = c(length(sel_cedrela$vars_cor),
                                length(sel_handroanthus$vars_cor)),
  N_vars_pos_VIF            = c(length(sel_cedrela$vars_final),
                                length(sel_handroanthus$vars_final)),
  AUC_treino_media          = c(round(resultado_cedrela$auc_treino["media"], 3),
                                round(resultado_handroanthus$auc_treino["media"], 3)),
  AUC_treino_sd             = c(round(resultado_cedrela$auc_treino["sd"], 3),
                                round(resultado_handroanthus$auc_treino["sd"], 3)),
  AUC_teste_media           = c(round(resultado_cedrela$auc_teste["media"], 3),
                                round(resultado_handroanthus$auc_teste["media"], 3)),
  AUC_teste_sd              = c(round(resultado_cedrela$auc_teste["sd"], 3),
                                round(resultado_handroanthus$auc_teste["sd"], 3)),
  AUC_geral                 = c(round(resultado_cedrela$auc_geral, 3),
                                round(resultado_handroanthus$auc_geral, 3)),
  Limiar_maxSSS             = c(round(tss_cedrela$limiar, 3),
                                round(tss_handroanthus$limiar, 3)),
  TSS_no_maxSSS             = c(round(tss_cedrela$TSS, 3),
                                round(tss_handroanthus$TSS, 3)),
  Sensibilidade             = c(round(tss_cedrela$sensibilidade, 3),
                                round(tss_handroanthus$sensibilidade, 3)),
  Especificidade            = c(round(tss_cedrela$especificidade, 3),
                                round(tss_handroanthus$especificidade, 3)),
  Omissao                   = c(round(tss_cedrela$omissao, 3),
                                round(tss_handroanthus$omissao, 3)),
  Area_adequada_km2_binario = c(
    sum(areasbin_cedrela$area_km2[areasbin_cedrela$valor == 1]),
    sum(areasbin_handroanthus$area_km2[areasbin_handroanthus$valor == 1])
  ),
  Pct_adequada_binario      = c(
    areasbin_cedrela$pct[areasbin_cedrela$valor == 1],
    areasbin_handroanthus$pct[areasbin_handroanthus$valor == 1]
  ),
  Area_medio_alta_km2       = c(
    sum(areas4_cedrela$area_km2[areas4_cedrela$valor %in% c(2, 3)]),
    sum(areas4_handroanthus$area_km2[areas4_handroanthus$valor %in% c(2, 3)])
  ),
  Pct_medio_alta            = c(
    sum(areas4_cedrela$pct[areas4_cedrela$valor %in% c(2, 3)]),
    sum(areas4_handroanthus$pct[areas4_handroanthus$valor %in% c(2, 3)])
  ),
  Vars_finais = c(paste(sel_cedrela$vars_final,      collapse = "; "),
                  paste(sel_handroanthus$vars_final, collapse = "; "))
)

print(resumo)

writexl::write_xlsx(
  list(
    resumo                = resumo,
    areas4_cedrela        = areas4_cedrela,
    areas4_handroanthus   = areas4_handroanthus,
    areasbin_cedrela      = areasbin_cedrela,
    areasbin_handroanthus = areasbin_handroanthus,
    vars_cedrela          = data.frame(variavel = sel_cedrela$vars_final),
    vars_handroanthus     = data.frame(variavel = sel_handroanthus$vars_final)
  ),
  file.path(out_dir, "resumo_maxent.xlsx")
)


#'
#' # 15. Referencia rapida: valores-alvo do artigo (apenas para H. serratifolius
#'      em escala continental). Aqui os valores serao diferentes porque:
#'      (a) restringimos ao Brasil; (b) tambem modelamos Cedrela.
#'
#'   AUC treino (artigo): 0.826 (SD 0.007)
#'   AUC teste  (artigo): 0.798 (SD 0.014)
#'   AUC geral  (artigo): 0.812
#'   TSS no maxSSS (artigo): 0.436 (SD 0.022)
#'   Limiar maxSSS (artigo): 0.53
#'   Area binaria adequada (artigo): 4.009.617 km^2 (~28%)
#'   Area medio-alta (artigo): ~3.976.801 km^2 (~27%)
#'
cat("\n=========================================================\n")
cat("Modelagem concluida. Saidas em: ", out_dir, "\n", sep = "")
cat("Principais arquivos:\n")
cat("  - *_adequabilidade_continua.tif    (mapa Cloglog 0-1)\n")
cat("  - *_4_classes.tif                  (estilo Fig. 4 do artigo)\n")
cat("  - *_binario_maxSSS.tif             (estilo Fig. 5 do artigo)\n")
cat("  - modelo_*.rds, tss_*.rds          (objetos para reuso)\n")
cat("  - resumo_maxent.xlsx               (Tab. comparavel ao artigo)\n")
cat("=========================================================\n")

save.image("data/resultados.RData")
#load("data/resultados.RData")