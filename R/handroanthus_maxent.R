# =========================================================================
# Modelos de nicho ecológico para Handroanthus no Brasil - MaxEnt
# =========================================================================

# =========================================================================
# 1. CONFIGURAÇÃO INICIAL
# =========================================================================

# Limpar ambiente, se desejar
# rm(list = ls())
# gc(reset = TRUE)
# graphics.off()

# Aumentar memória disponível para Java/MaxEnt
options(java.parameters = "-Xmx8g")

# Instalar e carregar pacotes
if(!require(readxl))    install.packages("readxl",    dependencies = TRUE)
if(!require(writexl))   install.packages("writexl",   dependencies = TRUE)
if(!require(tidyverse)) install.packages("tidyverse", dependencies = TRUE)
if(!require(terra))     install.packages("terra",     dependencies = TRUE)
if(!require(sf))        install.packages("sf",        dependencies = TRUE)
if(!require(sp))        install.packages("sp",        dependencies = TRUE)
if(!require(here))      install.packages("here",      dependencies = TRUE)
if(!require(dismo))     install.packages("dismo",     dependencies = TRUE)
if(!require(raster))    install.packages("raster",    dependencies = TRUE)
if(!require(rJava))     install.packages("rJava",     dependencies = TRUE)
if(!require(maps))      install.packages("maps",      dependencies = TRUE)
if(!require(geobr))     install.packages("geobr",     dependencies = TRUE)
if(!require(usdm))      install.packages("usdm",      dependencies = TRUE)
if(!require(spThin))    install.packages("spThin",    dependencies = TRUE)

library(readxl)
library(writexl)
library(tidyverse)
library(terra)
library(sf)
library(sp)
library(here)
library(dismo)
library(raster)
library(rJava)
library(maps)
library(geobr)
library(usdm)
library(spThin)
library(ggplot2)

# Semente para reprodutibilidade
set.seed(1350)

# Verificar se o arquivo maxent.jar está disponível
maxent_jar <- file.path(system.file("java", package = "dismo"), "maxent.jar")

if(!file.exists(maxent_jar)){
  warning("maxent.jar não foi encontrado em: ", maxent_jar,
          "\nBaixe o arquivo MaxEnt e coloque na pasta java do pacote dismo.")
}

# Diretórios de saída
out_dir <- here::here("maxent_output")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


out_handro  <- file.path(out_dir, "maxent_handroanthus")
if(!dir.exists(out_handro))  dir.create(out_handro,  recursive = TRUE)



# =========================================================================
# 2. CARREGAMENTO DOS DADOS E RASTERS
# =========================================================================

# Dados de ocorrência + variáveis ambientais já associadas
handroanthus <- read.csv(here::here("data", "handroanthus_var_amb.csv"))
head(handroanthus)
str(handroanthus)

# Carregar rasters ambientais
# Mantive o padrão do seu script original.
raster_files <- list.files(
  here::here("rasters"),
  #pattern = "32\\.vapor|33\\.solar\\.tif$",
  pattern = "\\.tif$",
  full.names = TRUE
)



# Excluir arquivos consolidados ou rasters regionais indesejados
excluir_pattern <- "Brazil_env_stack|Amazon_stack|Caatinga_elev_corridor|elev_corridor_minus1|48.density"
raster_files <- raster_files[!grepl(excluir_pattern, raster_files)]

# Stack em terra, útil para inspeção
r_stack <- terra::rast(raster_files)
names(r_stack) <- tools::file_path_sans_ext(basename(raster_files))

# Stack em raster, usado pelo dismo::maxent
r_stack_raster <- raster::stack(raster_files)
names(r_stack_raster) <- tools::file_path_sans_ext(basename(raster_files))

# Ver nomes originais
names(r_stack_raster)

# Remover prefixo do tipo X1., X10., X54., etc.
novos_nomes <- names(r_stack_raster)

novos_nomes <- gsub("^X[0-9]+\\.", "", novos_nomes)

# Aplicar nomes no RasterStack
names(r_stack_raster) <- novos_nomes

# Aplicar também no SpatRaster, se estiver usando terra
names(r_stack) <- novos_nomes

# Conferir
names(r_stack_raster)

cat("Número de camadas ambientais:", raster::nlayers(r_stack_raster), "\n")
print(names(r_stack_raster))

# Extensão espacial do estudo
ext_brasil <- raster::extent(r_stack_raster)
print(ext_brasil)

# Contorno do Brasil/Estados para mapas
brasil_uf <- tryCatch(
  geobr::read_state(year = 2020, showProgress = FALSE),
  error = function(e) NULL
)

if(!is.null(brasil_uf)){
  brasil_sp <- sf::as_Spatial(sf::st_geometry(brasil_uf))
} else {
  brasil_sp <- NULL
}


# =========================================================================
# 3. ESPÉCIE: handroanthus
# =========================================================================


# -------------------------------------------------------------------------
# 3.1 Preparação das ocorrências de handroanthus
# -------------------------------------------------------------------------

# Padronizar coordenadas
# Se o seu arquivo tem colunas x e y, usamos x = longitude e y = latitude.
handroanthus_pontos <- data.frame(
  lon = handroanthus$x,
  lat = handroanthus$y
)

# Remover valores ausentes
handroanthus_pontos <- handroanthus_pontos[complete.cases(handroanthus_pontos), ]

# Remover coordenadas duplicadas
handroanthus_pontos <- handroanthus_pontos[!duplicated(handroanthus_pontos), ]

# Remover coordenadas fora dos limites plausíveis
handroanthus_pontos <- handroanthus_pontos[
  handroanthus_pontos$lon >= -180 & handroanthus_pontos$lon <= 180 &
    handroanthus_pontos$lat >= -90 & handroanthus_pontos$lat <= 90,
]

cat("handroanthus - pontos antes do thinning:", nrow(handroanthus_pontos), "\n")

# Preparar tabela para spThin
handroanthus_thin_input <- handroanthus_pontos
handroanthus_thin_input$especie <- "handroanthus"

# Filtro espacial de 5 km
set.seed(1350)

handroanthus_thin_list <- spThin::thin(
  loc.data = handroanthus_thin_input,
  lat.col = "lat",
  long.col = "lon",
  spec.col = "especie",
  thin.par = 5,
  reps = 50,
  locs.thinned.list.return = TRUE,
  write.files = FALSE,
  write.log.file = FALSE,
  verbose = FALSE
)

# Selecionar a repetição que reteve o maior número de pontos
handroanthus_n_por_rep <- sapply(handroanthus_thin_list, nrow)
handroanthus_melhor_rep <- which.max(handroanthus_n_por_rep)
handroanthus_pts <- handroanthus_thin_list[[handroanthus_melhor_rep]]
names(handroanthus_pts) <- c("lon", "lat")

cat("handroanthus - pontos após thinning de 5 km:", nrow(handroanthus_pts), "\n")

# Visualizar pontos
plot(
  handroanthus_pts$lon,
  handroanthus_pts$lat,
  col = "darkgreen",
  pch = 19,
  cex = 0.5,
  xlab = "Longitude",
  ylab = "Latitude",
  main = paste0("handroanthus após thinning, n = ", nrow(handroanthus_pts))
)
maps::map(add = TRUE)


# -------------------------------------------------------------------------
# 3.2 Background e extração de valores ambientais para handroanthus
# -------------------------------------------------------------------------

set.seed(1350)

handroanthus_bg <- dismo::randomPoints(
  mask = r_stack_raster,
  n = 10000,
  ext = ext_brasil,
  extf = 1.1
)

colnames(handroanthus_bg) <- c("lon", "lat")
handroanthus_bg <- as.data.frame(handroanthus_bg)

# Extrair valores ambientais nos pontos de presença
handroanthus_vals_pres <- raster::extract(
  r_stack_raster,
  handroanthus_pts[, c("lon", "lat")],
  method = "simple"
)

# Extrair valores ambientais nos pontos de background
handroanthus_vals_bg <- raster::extract(
  r_stack_raster,
  handroanthus_bg[, c("lon", "lat")],
  method = "simple"
)

# Combinar presença + background para diagnóstico de colinearidade
handroanthus_vals_amb <- as.data.frame(rbind(handroanthus_vals_pres, handroanthus_vals_bg))
handroanthus_vals_amb <- handroanthus_vals_amb[complete.cases(handroanthus_vals_amb), ]

cat("handroanthus - linhas usadas para correlação/VIF:", nrow(handroanthus_vals_amb), "\n")


# -------------------------------------------------------------------------
# 3.3 Seleção de variáveis para handroanthus: correlação + VIF
# -------------------------------------------------------------------------

# Etapa 1: remover variáveis correlacionadas
# Observação: no texto metodológico, o ideal é descrever como Pearson |r| < 0.70.
# No código abaixo, usamos method = "pearson" para compatibilidade direta com isso.

handroanthus_sel_cor <- usdm::vifcor(
  handroanthus_vals_amb,
  th = 0.70,
  method = "pearson"
)

handroanthus_vars_cor <- as.character(handroanthus_sel_cor@results$Variables)

cat("handroanthus - número de variáveis após correlação:", length(handroanthus_vars_cor), "\n")
print(handroanthus_vars_cor)

# Etapa 2: remover variáveis com VIF alto
# No artigo citado, o limite descrito é VIF < 10.
# Use th = 10 para reproduzir esse critério.

handroanthus_sel_vif <- usdm::vifstep(
  handroanthus_vals_amb[, handroanthus_vars_cor],
  th = 10
)

handroanthus_vars_final <- as.character(handroanthus_sel_vif@results$Variables)

cat("handroanthus - número de variáveis após VIF:", length(handroanthus_vars_final), "\n")
print(handroanthus_vars_final)

# Stack final para handroanthus
stack_handroanthus_final <- raster::subset(r_stack_raster, handroanthus_vars_final)

cat("handroanthus - variáveis finais no RasterStack:\n")
print(names(stack_handroanthus_final))

# Extrair valores ambientais nos pontos de presença
handroanthus_vals_pres_final <- raster::extract(
  stack_handroanthus_final,
  handroanthus_pts[, c("lon", "lat")]
)

# Identificar pontos com pelo menos uma variável NA
handroanthus_pontos_com_NA <- !complete.cases(handroanthus_vals_pres_final)

# Quantidade e porcentagem
sum(handroanthus_pontos_com_NA)
100 * mean(handroanthus_pontos_com_NA)

# Pontos problemáticos
handroanthus_pts_NA <- handroanthus_pts[handroanthus_pontos_com_NA, ]

head(handroanthus_pts_NA)


# Número de NAs por variável ambiental
handroanthus_na_por_variavel <- colSums(is.na(handroanthus_vals_pres_final))

# Mostrar somente variáveis com NA
handroanthus_na_por_variavel <- handroanthus_na_por_variavel[handroanthus_na_por_variavel > 0]

# Ordenar da mais problemática para a menos problemática
sort(handroanthus_na_por_variavel, decreasing = TRUE)


plot(
  stack_handroanthus_final[[1]],
  main = "Pontos com e sem NA nas variáveis ambientais"
)

points(
  handroanthus_pts$lon,
  handroanthus_pts$lat,
  pch = 20,
  col = "blue",
  cex = 0.5
)

points(
  handroanthus_pts_NA$lon,
  handroanthus_pts_NA$lat,
  pch = 20,
  col = "red",
  cex = 0.8
)

legend(
  "bottomright",
  legend = c("Pontos válidos", "Pontos com NA"),
  col = c("blue", "red"),
  pch = 20,
  bty = "n"
)

# Manter apenas pontos sem NA em todas as variáveis ambientais
handroanthus_pts_modelo <- handroanthus_pts[!handroanthus_pontos_com_NA, ]

cat("handroanthus - pontos após thinning:", nrow(handroanthus_pts), "\n")
cat("handroanthus - pontos usados no MaxEnt:", nrow(handroanthus_pts_modelo), "\n")
cat("handroanthus - pontos removidos por NA:", sum(handroanthus_pontos_com_NA), "\n")

# -------------------------------------------------------------------------
# 3.4 Ajuste MaxEnt para handroanthus
# -------------------------------------------------------------------------

# Configuração final recomendada:
# - validação cruzada com 10 folds
# - 500 iterações máximas
# - saída cloglog
# - jackknife
# - curvas de resposta
# - regra Maximum training sensitivity plus specificity

handroanthus_args_maxent <- c(
  "replicates=10",
  "maximumiterations=500",
  "replicatetype=crossvalidate",
  "outputformat=cloglog",
  "betamultiplier=1",
  "threshold=true",
  "jackknife=true",
  "responsecurves=true",
  "writeplotdata=true",
  "writebackgroundpredictions=false",
  "applythresholdrule=Maximum training sensitivity plus specificity"
)

# Rodar modelo

cat("\n[6/8] Ajustando MaxEnt para handroanthus...\n")
cat("Início:", format(Sys.time(), "%d/%m/%Y %H:%M:%S"), "\n")
tempo_handroanthus <- system.time({
  
  set.seed(1350)
  handroanthus_modelo <- dismo::maxent(
    x = stack_handroanthus_final,
    p = handroanthus_pts_modelo[, c("lon", "lat")],
    args = handroanthus_args_maxent,
    path = out_handro,
    silent = FALSE
  )
  
})

cat("Fim:", format(Sys.time(), "%d/%m/%Y %H:%M:%S"), "\n")
cat("Tempo total do ajuste MaxEnt - handroanthus:\n")
print(tempo_handroanthus)

print(handroanthus_modelo)

#--------------------------------------------------------------------------------------
# Gráfico de importancia de variáveis
# -------------------------------------------------------------------------------------

# Extrair resultados do modelo
handroanthus_resultados <- handroanthus_modelo@results

# Linhas com percentual de contribuição
linhas_contrib <- grep("\\.contribution$", rownames(handroanthus_resultados), value = TRUE)

# Nomes das variáveis
vars <- gsub("\\.contribution$", "", linhas_contrib)

# Média da contribuição entre réplicas
contrib_media <- rowMeans(
  as.matrix(handroanthus_resultados[linhas_contrib, , drop = FALSE]),
  na.rm = TRUE
)

# Linhas com importância por permutação
linhas_perm <- paste0(vars, ".permutation.importance")

perm_media <- rowMeans(
  as.matrix(handroanthus_resultados[linhas_perm, , drop = FALSE]),
  na.rm = TRUE
)

# Montar tabela
imp_handroanthus <- data.frame(
  variavel = vars,
  contribuicao = contrib_media,
  importancia_permutacao = perm_media
)

# Ordenar
imp_handroanthus <- imp_handroanthus[order(imp_handroanthus$contribuicao, decreasing = TRUE), ]

imp_handroanthus




ggplot(imp_handroanthus,
       aes(x = reorder(variavel, contribuicao),
           y = contribuicao)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Variável ambiental",
    y = "Contribuição (%)",
    title = "Importância das variáveis - handroanthus",
    subtitle = "Percentual de contribuição no modelo MaxEnt"
  ) +
  theme_bw()



library(tidyr)
library(ggplot2)

imp_handroanthus_long <- imp_handroanthus |>
  pivot_longer(
    cols = c(contribuicao, importancia_permutacao),
    names_to = "metrica",
    values_to = "valor"
  )

ggplot(imp_handroanthus_long,
       aes(x = reorder(variavel, valor),
           y = valor,
           fill = metrica)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    x = "Variável ambiental",
    y = "Importância (%)",
    fill = "Métrica",
    title = "Importância das variáveis - handroanthus"
  ) +
  theme_bw()


# -------------------------------------------------------------------------
# 3.5 AUC com validação cruzada para handroanthus
# -------------------------------------------------------------------------

# Verificar estrutura
rownames(handroanthus_resultados)
colnames(handroanthus_resultados)

# -------------------------------------------------------------------------
# Identificar colunas das réplicas da validação cruzada
# -------------------------------------------------------------------------

# Em geral, as colunas aparecem como:
# species_0, species_1, ..., species_9
# ou handroanthus_0, handroanthus_1, ..., handroanthus_9

handroanthus_cols <- colnames(handroanthus_resultados)

handroanthus_cols_rep <- handroanthus_cols[
  grepl("_[0-9]+$", handroanthus_cols)
]

# Caso o padrão acima não funcione, remover colunas de média/desvio
if(length(handroanthus_cols_rep) == 0){
  handroanthus_cols_rep <- handroanthus_cols[
    !grepl("average|avg|std|stddev|sd", handroanthus_cols, ignore.case = TRUE)
  ]
}

cat("Colunas das réplicas usadas no cálculo da AUC:\n")
print(handroanthus_cols_rep)

# -------------------------------------------------------------------------
# Identificar linhas exatas de AUC
# -------------------------------------------------------------------------

linha_auc_treino <- which(rownames(handroanthus_resultados) == "Training.AUC")
linha_auc_teste  <- which(rownames(handroanthus_resultados) == "Test.AUC")

if(length(linha_auc_treino) == 0){
  stop("Linha 'Training.AUC' não encontrada em handroanthus_modelo@results.")
}

if(length(linha_auc_teste) == 0){
  stop("Linha 'Test.AUC' não encontrada em handroanthus_modelo@results.")
}

# -------------------------------------------------------------------------
# Extrair AUC por fold
# -------------------------------------------------------------------------

handroanthus_auc_cv <- data.frame(
  fold = seq_along(handroanthus_cols_rep),
  replica = handroanthus_cols_rep,
  AUC_treino = as.numeric(handroanthus_resultados[linha_auc_treino, handroanthus_cols_rep]),
  AUC_teste  = as.numeric(handroanthus_resultados[linha_auc_teste,  handroanthus_cols_rep])
)

# Diferença treino - teste como diagnóstico de possível sobreajuste
handroanthus_auc_cv$delta_treino_teste <- handroanthus_auc_cv$AUC_treino - handroanthus_auc_cv$AUC_teste

print(handroanthus_auc_cv)

# -------------------------------------------------------------------------
# Resumo da AUC por validação cruzada
# -------------------------------------------------------------------------

handroanthus_auc_resumo <- data.frame(
  metrica = c("AUC treino", "AUC teste - validação cruzada", "Diferença treino - teste"),
  media = c(
    mean(handroanthus_auc_cv$AUC_treino, na.rm = TRUE),
    mean(handroanthus_auc_cv$AUC_teste, na.rm = TRUE),
    mean(handroanthus_auc_cv$delta_treino_teste, na.rm = TRUE)
  ),
  desvio_padrao = c(
    sd(handroanthus_auc_cv$AUC_treino, na.rm = TRUE),
    sd(handroanthus_auc_cv$AUC_teste, na.rm = TRUE),
    sd(handroanthus_auc_cv$delta_treino_teste, na.rm = TRUE)
  ),
  minimo = c(
    min(handroanthus_auc_cv$AUC_treino, na.rm = TRUE),
    min(handroanthus_auc_cv$AUC_teste, na.rm = TRUE),
    min(handroanthus_auc_cv$delta_treino_teste, na.rm = TRUE)
  ),
  maximo = c(
    max(handroanthus_auc_cv$AUC_treino, na.rm = TRUE),
    max(handroanthus_auc_cv$AUC_teste, na.rm = TRUE),
    max(handroanthus_auc_cv$delta_treino_teste, na.rm = TRUE)
  ),
  n_folds = c(
    sum(!is.na(handroanthus_auc_cv$AUC_treino)),
    sum(!is.na(handroanthus_auc_cv$AUC_teste)),
    sum(!is.na(handroanthus_auc_cv$delta_treino_teste))
  )
)

print(handroanthus_auc_resumo)

# -------------------------------------------------------------------------
# AUC principal a reportar no artigo
# -------------------------------------------------------------------------

handroanthus_auc_treino_media <- mean(handroanthus_auc_cv$AUC_treino, na.rm = TRUE)
handroanthus_auc_treino_sd    <- sd(handroanthus_auc_cv$AUC_treino, na.rm = TRUE)

handroanthus_auc_teste_media <- mean(handroanthus_auc_cv$AUC_teste, na.rm = TRUE)
handroanthus_auc_teste_sd    <- sd(handroanthus_auc_cv$AUC_teste, na.rm = TRUE)

handroanthus_delta_auc_media <- mean(handroanthus_auc_cv$delta_treino_teste, na.rm = TRUE)

cat("\nhandroanthus - AUC treino:",
    round(handroanthus_auc_treino_media, 3),
    "±",
    round(handroanthus_auc_treino_sd, 3),
    "\n")

cat("handroanthus - AUC teste por validação cruzada:",
    round(handroanthus_auc_teste_media, 3),
    "±",
    round(handroanthus_auc_teste_sd, 3),
    "\n")

cat("handroanthus - diferença média AUC treino - teste:",
    round(handroanthus_delta_auc_media, 3),
    "\n")

# -------------------------------------------------------------------------
# 3.6 Predição espacial contínua para handroanthus
# -------------------------------------------------------------------------

cat("\nGerando predição espacial para handroanthus...\n")

handroanthus_pred_bruto <- predict(
  handroanthus_modelo,
  stack_handroanthus_final,
  progress = "text"
)

class(handroanthus_pred_bruto)

# Se o MaxEnt retornar uma camada por fold,
# calcular o mapa médio e o desvio-padrão entre folds
if(inherits(handroanthus_pred_bruto, "RasterStack") | 
   inherits(handroanthus_pred_bruto, "RasterBrick")){
  
  handroanthus_pred <- raster::calc(
    handroanthus_pred_bruto,
    fun = mean,
    na.rm = TRUE
  )
  
  handroanthus_pred_sd <- raster::calc(
    handroanthus_pred_bruto,
    fun = sd,
    na.rm = TRUE
  )
  
} else {
  
  handroanthus_pred <- handroanthus_pred_bruto
  handroanthus_pred_sd <- NULL
  
}

# -------------------------------------------------------------------------
# Mapa contínuo médio em ggplot2
# -------------------------------------------------------------------------

library(ggplot2)
library(sf)

# Converter raster para data.frame
handroanthus_pred_df <- raster::as.data.frame(
  handroanthus_pred,
  xy = TRUE,
  na.rm = TRUE
)

# Padronizar nome da coluna de predição
names(handroanthus_pred_df)[3] <- "adequabilidade"

# Converter contorno do Brasil para sf, se disponível
brasil_sf <- NULL

if(exists("brasil_uf") && !is.null(brasil_uf)){
  brasil_sf <- sf::st_as_sf(brasil_uf)
}

if(is.null(brasil_sf) && exists("brasil_sp") && !is.null(brasil_sp)){
  brasil_sf <- sf::st_as_sf(brasil_sp)
}

# Extensão do raster
ext_handroanthus <- raster::extent(handroanthus_pred)

# Camada do Brasil
camada_brasil <- if(!is.null(brasil_sf)){
  geom_sf(
    data = brasil_sf,
    fill = NA,
    color = "gray30",
    linewidth = 0.3,
    inherit.aes = FALSE
  )
} else {
  NULL
}

# Mapa em ggplot2
mapa_handroanthus_gg <- ggplot() +
  geom_raster(
    data = handroanthus_pred_df,
    aes(
      x = x,
      y = y,
      fill = adequabilidade
    )
  ) +
  camada_brasil +
  geom_point(
    data = handroanthus_pts_modelo,
    aes(
      x = lon,
      y = lat
    ),
    shape = 21,
    fill = "red",
    color = "black",
    size = 1.2,
    stroke = 0.25
  ) +
  scale_fill_viridis_c(
    name = "Adequabilidade",
    limits = c(0, 1),
    na.value = "transparent"
  ) +
  coord_sf(
    xlim = c(ext_handroanthus@xmin, ext_handroanthus@xmax),
    ylim = c(ext_handroanthus@ymin, ext_handroanthus@ymax),
    expand = FALSE
  ) +
  labs(
    title = expression(italic("handroanthus") ~ "- adequabilidade MaxEnt"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid = element_line(linewidth = 0.2, color = "gray85")
  )

mapa_handroanthus_gg

#Salvar figura
ggsave(
  filename =  "figuras/handroanthus_adequabilidade.png",
  plot = mapa_handroanthus_gg,
  width = 8,
  height = 6,
  dpi = 300
)


# -------------------------------------------------------------------------
# Mapa de desvio-padrão entre folds em ggplot2
# -------------------------------------------------------------------------

if(!is.null(handroanthus_pred_sd)){
  
  # Converter raster para data.frame
  handroanthus_pred_sd_df <- raster::as.data.frame(
    handroanthus_pred_sd,
    xy = TRUE,
    na.rm = TRUE
  )
  
  # Padronizar nome da coluna
  names(handroanthus_pred_sd_df)[3] <- "desvio_padrao"
  
  # Converter contorno do Brasil para sf, se disponível
  brasil_sf <- NULL
  
  if(exists("brasil_uf") && !is.null(brasil_uf)){
    brasil_sf <- sf::st_as_sf(brasil_uf)
  }
  
  if(is.null(brasil_sf) && exists("brasil_sp") && !is.null(brasil_sp)){
    brasil_sf <- sf::st_as_sf(brasil_sp)
  }
  
  # Extensão do raster
  ext_handroanthus_sd <- raster::extent(handroanthus_pred_sd)
  
  # Camada do Brasil
  camada_brasil_sd <- if(!is.null(brasil_sf)){
    geom_sf(
      data = brasil_sf,
      fill = NA,
      color = "gray30",
      linewidth = 0.3,
      inherit.aes = FALSE
    )
  } else {
    NULL
  }
  
  # Mapa em ggplot2
  mapa_handroanthus_sd_gg <- ggplot() +
    geom_raster(
      data = handroanthus_pred_sd_df,
      aes(
        x = x,
        y = y,
        fill = desvio_padrao
      )
    ) +
    camada_brasil_sd +
    scale_fill_viridis_c(
      name = "Desvio-padrão\nentre folds",
      na.value = "transparent"
    ) +
    coord_sf(
      xlim = c(ext_handroanthus_sd@xmin, ext_handroanthus_sd@xmax),
      ylim = c(ext_handroanthus_sd@ymin, ext_handroanthus_sd@ymax),
      expand = FALSE
    ) +
    labs(
      title = expression(italic("handroanthus") ~ "- desvio-padrão entre folds"),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      panel.grid = element_line(linewidth = 0.2, color = "gray85")
    )
  
  mapa_handroanthus_sd_gg
}


if(!is.null(handroanthus_pred_sd)){
  
  ggsave(
    filename = "figuras/handroanthus_desvio_padrao.png",
    plot = mapa_handroanthus_sd_gg,
    width = 8,
    height = 6,
    dpi = 300
  )
  
}

# -------------------------------------------------------------------------
# 3.7 Limiar maxSSS e TSS para handroanthus
# -------------------------------------------------------------------------

# Verificar métricas de threshold disponíveis
grep(
  "Maximum.training.sensitivity.plus.specificity",
  rownames(handroanthus_resultados),
  value = TRUE
)

# Buscar preferencialmente o limiar em Cloglog
handroanthus_linhas_limiar <- grep(
  "Maximum.training.sensitivity.plus.specificity.Cloglog.threshold",
  rownames(handroanthus_resultados)
)

# Caso o nome exato não exista, buscar qualquer threshold maxSSS
if(length(handroanthus_linhas_limiar) == 0){
  handroanthus_linhas_limiar <- grep(
    "Maximum.training.sensitivity.plus.specificity.*threshold",
    rownames(handroanthus_resultados)
  )
}

if(length(handroanthus_linhas_limiar) == 0){
  stop("Limiar Maximum training sensitivity plus specificity não encontrado.")
}

# Extrair limiares apenas das réplicas/folds
handroanthus_limiar_vals <- as.numeric(
  handroanthus_resultados[handroanthus_linhas_limiar, handroanthus_cols_rep, drop = FALSE]
)

handroanthus_limiar_vals <- handroanthus_limiar_vals[!is.na(handroanthus_limiar_vals)]

handroanthus_limiar_maxSSS <- mean(handroanthus_limiar_vals, na.rm = TRUE)
handroanthus_limiar_maxSSS_sd <- sd(handroanthus_limiar_vals, na.rm = TRUE)

cat("handroanthus - limiar maxSSS médio:",
    round(handroanthus_limiar_maxSSS, 3),
    "±",
    round(handroanthus_limiar_maxSSS_sd, 3),
    "\n")

# -------------------------------------------------------------------------
# Avaliar sensibilidade, especificidade, TSS e omissão no limiar maxSSS
# -------------------------------------------------------------------------

set.seed(1350)

handroanthus_bg_tss <- dismo::randomPoints(
  mask = stack_handroanthus_final,
  n = 10000,
  ext = ext_brasil
)

colnames(handroanthus_bg_tss) <- c("lon", "lat")
handroanthus_bg_tss <- as.data.frame(handroanthus_bg_tss)

# Extrair adequabilidade nos pontos de presença usados no modelo
handroanthus_s_pres <- raster::extract(
  handroanthus_pred,
  handroanthus_pts_modelo[, c("lon", "lat")]
)

# Extrair adequabilidade nos pontos de background
handroanthus_s_bg <- raster::extract(
  handroanthus_pred,
  handroanthus_bg_tss[, c("lon", "lat")]
)

# Remover NAs
handroanthus_s_pres <- handroanthus_s_pres[!is.na(handroanthus_s_pres)]
handroanthus_s_bg <- handroanthus_s_bg[!is.na(handroanthus_s_bg)]

# Matriz de confusão usando presença e background
handroanthus_TP <- sum(handroanthus_s_pres >= handroanthus_limiar_maxSSS)
handroanthus_FN <- sum(handroanthus_s_pres <  handroanthus_limiar_maxSSS)

handroanthus_FP <- sum(handroanthus_s_bg >= handroanthus_limiar_maxSSS)
handroanthus_TN <- sum(handroanthus_s_bg <  handroanthus_limiar_maxSSS)

# Métricas
handroanthus_sensibilidade <- handroanthus_TP / (handroanthus_TP + handroanthus_FN)
handroanthus_especificidade <- handroanthus_TN / (handroanthus_TN + handroanthus_FP)

handroanthus_TSS <- handroanthus_sensibilidade + handroanthus_especificidade - 1
handroanthus_omissao <- handroanthus_FN / (handroanthus_TP + handroanthus_FN)

handroanthus_metricas_threshold <- data.frame(
  especie = "handroanthus",
  limiar_maxSSS = handroanthus_limiar_maxSSS,
  sensibilidade = handroanthus_sensibilidade,
  especificidade = handroanthus_especificidade,
  TSS = handroanthus_TSS,
  omissao = handroanthus_omissao,
  TP = handroanthus_TP,
  FN = handroanthus_FN,
  FP = handroanthus_FP,
  TN = handroanthus_TN
)

print(handroanthus_metricas_threshold)

cat("handroanthus - sensibilidade:", round(handroanthus_sensibilidade, 3), "\n")
cat("handroanthus - especificidade:", round(handroanthus_especificidade, 3), "\n")
cat("handroanthus - TSS:", round(handroanthus_TSS, 3), "\n")
cat("handroanthus - omissão:", round(handroanthus_omissao, 3), "\n")


# -------------------------------------------------------------------------
# 3.8 Mapas classificados para handroanthus em ggplot2
# -------------------------------------------------------------------------

library(ggplot2)
library(sf)
library(dplyr)

# -------------------------------------------------------------------------
# Reclassificação: 4 classes de adequabilidade
# -------------------------------------------------------------------------

handroanthus_rcl_4classes <- matrix(
  c(
    -Inf, 0.265, 0,
    0.265, 0.530, 1,
    0.530, 0.765, 2,
    0.765, Inf,   3
  ),
  ncol = 3,
  byrow = TRUE
)

handroanthus_mapa4 <- raster::reclassify(
  handroanthus_pred,
  rcl = handroanthus_rcl_4classes
)

# -------------------------------------------------------------------------
# Reclassificação: mapa binário maxSSS
# -------------------------------------------------------------------------

handroanthus_rcl_bin <- matrix(
  c(
    -Inf, handroanthus_limiar_maxSSS, 0,
    handroanthus_limiar_maxSSS, Inf, 1
  ),
  ncol = 3,
  byrow = TRUE
)

handroanthus_mapa_bin <- raster::reclassify(
  handroanthus_pred,
  rcl = handroanthus_rcl_bin
)

# -------------------------------------------------------------------------
# Preparar contorno do Brasil para ggplot2
# -------------------------------------------------------------------------

brasil_sf <- NULL

if(exists("brasil_uf") && !is.null(brasil_uf)){
  brasil_sf <- sf::st_as_sf(brasil_uf)
}

if(is.null(brasil_sf) && exists("brasil_sp") && !is.null(brasil_sp)){
  brasil_sf <- sf::st_as_sf(brasil_sp)
}

camada_brasil <- if(!is.null(brasil_sf)){
  geom_sf(
    data = brasil_sf,
    fill = NA,
    color = "gray30",
    linewidth = 0.3,
    inherit.aes = FALSE
  )
} else {
  NULL
}

# Extensão espacial
ext_handroanthus <- raster::extent(handroanthus_pred)

# -------------------------------------------------------------------------
# Mapa 4 classes em ggplot2
# -------------------------------------------------------------------------

handroanthus_mapa4_df <- raster::as.data.frame(
  handroanthus_mapa4,
  xy = TRUE,
  na.rm = TRUE
)

names(handroanthus_mapa4_df)[3] <- "classe_valor"

handroanthus_mapa4_df <- handroanthus_mapa4_df |>
  mutate(
    classe = factor(
      classe_valor,
      levels = c(0, 1, 2, 3),
      labels = c(
        "Inadequado (<0.265)",
        "Baixa (0.265-0.530)",
        "Média (0.530-0.765)",
        "Alta (>=0.765)"
      )
    )
  )

cores_4cl <- c(
  "Inadequado (<0.265)" = "#f0f0f0",
  "Baixa (0.265-0.530)" = "#fee08b",
  "Média (0.530-0.765)" = "#a6d96a",
  "Alta (>=0.765)" = "#1a9850"
)

mapa_handroanthus_4classes_gg <- ggplot() +
  geom_raster(
    data = handroanthus_mapa4_df,
    aes(
      x = x,
      y = y,
      fill = classe
    )
  ) +
  camada_brasil +
  scale_fill_manual(
    name = "Classe de\nadequabilidade",
    values = cores_4cl,
    drop = FALSE,
    na.value = "transparent"
  ) +
  coord_sf(
    xlim = c(ext_handroanthus@xmin, ext_handroanthus@xmax),
    ylim = c(ext_handroanthus@ymin, ext_handroanthus@ymax),
    expand = FALSE
  ) +
  labs(
    title = expression(italic("handroanthus") ~ "- classes de adequabilidade"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid = element_line(linewidth = 0.2, color = "gray85")
  )

mapa_handroanthus_4classes_gg

#Salvar
ggsave(
  filename = "figuras/handroanthus_4_classes_ggplot.png",
  plot = mapa_handroanthus_4classes_gg,
  width = 8,
  height = 6,
  dpi = 300
)


# -------------------------------------------------------------------------
# Mapa binário em ggplot2
# -------------------------------------------------------------------------

handroanthus_mapa_bin_df <- raster::as.data.frame(
  handroanthus_mapa_bin,
  xy = TRUE,
  na.rm = TRUE
)

names(handroanthus_mapa_bin_df)[3] <- "classe_valor"

handroanthus_mapa_bin_df <- handroanthus_mapa_bin_df |>
  mutate(
    classe = factor(
      classe_valor,
      levels = c(0, 1),
      labels = c("Inadequado", "Adequado")
    )
  )

cores_bin <- c(
  "Inadequado" = "#d9d9d9",
  "Adequado" = "#c2185b"
)

mapa_handroanthus_binario_gg <- ggplot() +
  geom_raster(
    data = handroanthus_mapa_bin_df,
    aes(
      x = x,
      y = y,
      fill = classe
    )
  ) +
  camada_brasil +
  geom_point(
    data = handroanthus_pts_modelo,
    aes(
      x = lon,
      y = lat
    ),
    shape = 21,
    fill = "yellow",
    color = "black",
    size = 1.1,
    stroke = 0.25
  ) +
  scale_fill_manual(
    name = "Classe",
    values = cores_bin,
    drop = FALSE,
    na.value = "transparent"
  ) +
  coord_sf(
    xlim = c(ext_handroanthus@xmin, ext_handroanthus@xmax),
    ylim = c(ext_handroanthus@ymin, ext_handroanthus@ymax),
    expand = FALSE
  ) +
  labs(
    title = paste0(
      "handroanthus - binário maxSSS = ",
      round(handroanthus_limiar_maxSSS, 3)
    ),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid = element_line(linewidth = 0.2, color = "gray85")
  )

mapa_handroanthus_binario_gg



#Salvar figura
ggsave(
  filename = "figuras/handroanthus_binario_maxSSS_ggplot.png",
  plot = mapa_handroanthus_binario_gg,
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------------
# 3.9 Área por classe para handroanthus
# -------------------------------------------------------------------------

# Área de cada célula em km²
handroanthus_area_pixel_4 <- raster::area(handroanthus_mapa4, na.rm = TRUE)

handroanthus_vals_classe4 <- raster::values(handroanthus_mapa4)
handroanthus_vals_area4 <- raster::values(handroanthus_area_pixel_4)

handroanthus_ok4 <- !is.na(handroanthus_vals_classe4) & !is.na(handroanthus_vals_area4)

handroanthus_vals_classe4 <- handroanthus_vals_classe4[handroanthus_ok4]
handroanthus_vals_area4 <- handroanthus_vals_area4[handroanthus_ok4]

handroanthus_area4_soma <- tapply(
  handroanthus_vals_area4,
  handroanthus_vals_classe4,
  sum
)

handroanthus_area4_total <- sum(handroanthus_area4_soma)

handroanthus_areas4 <- data.frame(
  valor = as.integer(names(handroanthus_area4_soma)),
  classe = c("Inadequado", "Baixa", "Média", "Alta")[
    as.integer(names(handroanthus_area4_soma)) + 1
  ],
  area_km2 = round(as.numeric(handroanthus_area4_soma), 0),
  pct = round(
    100 * as.numeric(handroanthus_area4_soma) / handroanthus_area4_total,
    2
  )
)

print(handroanthus_areas4)

# Área binária
handroanthus_area_pixel_bin <- raster::area(handroanthus_mapa_bin, na.rm = TRUE)

handroanthus_vals_bin <- raster::values(handroanthus_mapa_bin)
handroanthus_vals_area_bin <- raster::values(handroanthus_area_pixel_bin)

handroanthus_ok_bin <- !is.na(handroanthus_vals_bin) & !is.na(handroanthus_vals_area_bin)

handroanthus_vals_bin <- handroanthus_vals_bin[handroanthus_ok_bin]
handroanthus_vals_area_bin <- handroanthus_vals_area_bin[handroanthus_ok_bin]

handroanthus_area_bin_soma <- tapply(
  handroanthus_vals_area_bin,
  handroanthus_vals_bin,
  sum
)

handroanthus_area_bin_total <- sum(handroanthus_area_bin_soma)

handroanthus_areasbin <- data.frame(
  valor = as.integer(names(handroanthus_area_bin_soma)),
  classe = c("Inadequado", "Adequado")[
    as.integer(names(handroanthus_area_bin_soma)) + 1
  ],
  area_km2 = round(as.numeric(handroanthus_area_bin_soma), 0),
  pct = round(
    100 * as.numeric(handroanthus_area_bin_soma) / handroanthus_area_bin_total,
    2
  )
)

print(handroanthus_areasbin)


# -------------------------------------------------------------------------
# 3.10 Exportação dos resultados de handroanthus
# -------------------------------------------------------------------------

# Raster contínuo médio
raster::writeRaster(
  handroanthus_pred,
  file.path(out_dir, "handroanthus_adequabilidade_continua_media_cv.tif"),
  overwrite = TRUE
)

# Raster de desvio-padrão entre folds
if(!is.null(handroanthus_pred_sd)){
  
  raster::writeRaster(
    handroanthus_pred_sd,
    file.path(out_dir, "handroanthus_adequabilidade_sd_cv.tif"),
    overwrite = TRUE
  )
  
}

# Raster com 4 classes
raster::writeRaster(
  handroanthus_mapa4,
  file.path(out_dir, "handroanthus_4_classes.tif"),
  overwrite = TRUE,
  datatype = "INT1U"
)

# Raster binário maxSSS
raster::writeRaster(
  handroanthus_mapa_bin,
  file.path(out_dir, "handroanthus_binario_maxSSS.tif"),
  overwrite = TRUE,
  datatype = "INT1U"
)

# Tabelas de resultados
writexl::write_xlsx(
  list(
    AUC_por_fold = handroanthus_auc_cv,
    AUC_resumo = handroanthus_auc_resumo,
    Metricas_threshold = handroanthus_metricas_threshold,
    Area_4_classes = handroanthus_areas4,
    Area_binaria = handroanthus_areasbin
  ),
  path = file.path(out_dir, "handroanthus_resultados_maxent_cv.xlsx")
)

# Objetos R
saveRDS(
  handroanthus_modelo,
  file.path(out_dir, "modelo_handroanthus_cv.rds")
)

saveRDS(
  handroanthus_pred,
  file.path(out_dir, "predicao_handroanthus_media_cv.rds")
)

if(!is.null(handroanthus_pred_sd)){
  saveRDS(
    handroanthus_pred_sd,
    file.path(out_dir, "predicao_handroanthus_sd_cv.rds")
  )
}

# Salvar ambiente de trabalho
save.image("data/modelo_handroanthus_cv.RData")

