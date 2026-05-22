# =========================================================================
# Modelos de nicho ecológico para Cedrela e Handroanthus no Brasil - MaxEnt
# Versão simplificada, transparente e organizada por espécie
#
# Objetivo desta versão:
#   - evitar funções personalizadas;
#   - deixar o fluxo explícito;
#   - repetir os blocos por espécie;
#   - facilitar inspeção, depuração e adaptação manual.
#
# Estrutura geral:
#   1. Configuração inicial
#   2. Carregamento dos dados e rasters
#   3. Espécie 1: Cedrela
#   4. Espécie 2: Handroanthus
#   5. Tabelas e exportações finais
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

out_cedrela <- file.path(out_dir, "maxent_cedrela")
out_handro  <- file.path(out_dir, "maxent_handroanthus")

if(!dir.exists(out_cedrela)) dir.create(out_cedrela, recursive = TRUE)
if(!dir.exists(out_handro))  dir.create(out_handro,  recursive = TRUE)


# =========================================================================
# 2. CARREGAMENTO DOS DADOS E RASTERS
# =========================================================================

# Dados de ocorrência + variáveis ambientais já associadas
cedrela <- read.csv(here::here("data", "cedrela_br_var_amb.csv"))
handroanthus <- read.csv(here::here("data", "handroanthus_var_amb.csv"))

head(cedrela)
head(handroanthus)
str(cedrela)
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
# 3. ESPÉCIE 1: Cedrela
# =========================================================================


# -------------------------------------------------------------------------
# 3.1 Preparação das ocorrências de Cedrela
# -------------------------------------------------------------------------

# Padronizar coordenadas
# Se o seu arquivo tem colunas x e y, usamos x = longitude e y = latitude.
cedrela_pontos <- data.frame(
  lon = cedrela$x,
  lat = cedrela$y
)

# Remover valores ausentes
cedrela_pontos <- cedrela_pontos[complete.cases(cedrela_pontos), ]

# Remover coordenadas duplicadas
cedrela_pontos <- cedrela_pontos[!duplicated(cedrela_pontos), ]

# Remover coordenadas fora dos limites plausíveis
cedrela_pontos <- cedrela_pontos[
  cedrela_pontos$lon >= -180 & cedrela_pontos$lon <= 180 &
    cedrela_pontos$lat >= -90 & cedrela_pontos$lat <= 90,
]

cat("Cedrela - pontos antes do thinning:", nrow(cedrela_pontos), "\n")

# Preparar tabela para spThin
cedrela_thin_input <- cedrela_pontos
cedrela_thin_input$especie <- "Cedrela"

# Filtro espacial de 5 km
set.seed(1350)

cedrela_thin_list <- spThin::thin(
  loc.data = cedrela_thin_input,
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
cedrela_n_por_rep <- sapply(cedrela_thin_list, nrow)
cedrela_melhor_rep <- which.max(cedrela_n_por_rep)
cedrela_pts <- cedrela_thin_list[[cedrela_melhor_rep]]
names(cedrela_pts) <- c("lon", "lat")

cat("Cedrela - pontos após thinning de 5 km:", nrow(cedrela_pts), "\n")

# Visualizar pontos
plot(
  cedrela_pts$lon,
  cedrela_pts$lat,
  col = "darkgreen",
  pch = 19,
  cex = 0.5,
  xlab = "Longitude",
  ylab = "Latitude",
  main = paste0("Cedrela após thinning, n = ", nrow(cedrela_pts))
)
maps::map(add = TRUE)


# -------------------------------------------------------------------------
# 3.2 Background e extração de valores ambientais para Cedrela
# -------------------------------------------------------------------------

set.seed(1350)

cedrela_bg <- dismo::randomPoints(
  mask = r_stack_raster,
  n = 10000,
  ext = ext_brasil,
  extf = 1.1
)

colnames(cedrela_bg) <- c("lon", "lat")
cedrela_bg <- as.data.frame(cedrela_bg)

# Extrair valores ambientais nos pontos de presença
cedrela_vals_pres <- raster::extract(
  r_stack_raster,
  cedrela_pts[, c("lon", "lat")],
  method = "simple"
)

# Extrair valores ambientais nos pontos de background
cedrela_vals_bg <- raster::extract(
  r_stack_raster,
  cedrela_bg[, c("lon", "lat")],
  method = "simple"
)

# Combinar presença + background para diagnóstico de colinearidade
cedrela_vals_amb <- as.data.frame(rbind(cedrela_vals_pres, cedrela_vals_bg))
cedrela_vals_amb <- cedrela_vals_amb[complete.cases(cedrela_vals_amb), ]

cat("Cedrela - linhas usadas para correlação/VIF:", nrow(cedrela_vals_amb), "\n")


# -------------------------------------------------------------------------
# 3.3 Seleção de variáveis para Cedrela: correlação + VIF
# -------------------------------------------------------------------------

# Etapa 1: remover variáveis correlacionadas
# Observação: no texto metodológico, o ideal é descrever como Pearson |r| < 0.70.
# No código abaixo, usamos method = "pearson" para compatibilidade direta com isso.

cedrela_sel_cor <- usdm::vifcor(
  cedrela_vals_amb,
  th = 0.70,
  method = "pearson"
)

cedrela_vars_cor <- as.character(cedrela_sel_cor@results$Variables)

cat("Cedrela - número de variáveis após correlação:", length(cedrela_vars_cor), "\n")
print(cedrela_vars_cor)

# Etapa 2: remover variáveis com VIF alto
# No artigo citado, o limite descrito é VIF < 10.
# Use th = 10 para reproduzir esse critério.

cedrela_sel_vif <- usdm::vifstep(
  cedrela_vals_amb[, cedrela_vars_cor],
  th = 10
)

cedrela_vars_final <- as.character(cedrela_sel_vif@results$Variables)

cat("Cedrela - número de variáveis após VIF:", length(cedrela_vars_final), "\n")
print(cedrela_vars_final)

# Stack final para Cedrela
stack_cedrela_final <- raster::subset(r_stack_raster, cedrela_vars_final)

cat("Cedrela - variáveis finais no RasterStack:\n")
print(names(stack_cedrela_final))

# Extrair valores ambientais nos pontos de presença
cedrela_vals_pres_final <- raster::extract(
  stack_cedrela_final,
  cedrela_pts[, c("lon", "lat")]
)

# Identificar pontos com pelo menos uma variável NA
cedrela_pontos_com_NA <- !complete.cases(cedrela_vals_pres_final)

# Quantidade e porcentagem
sum(cedrela_pontos_com_NA)
100 * mean(cedrela_pontos_com_NA)

# Pontos problemáticos
cedrela_pts_NA <- cedrela_pts[cedrela_pontos_com_NA, ]

head(cedrela_pts_NA)


# Número de NAs por variável ambiental
cedrela_na_por_variavel <- colSums(is.na(cedrela_vals_pres_final))

# Mostrar somente variáveis com NA
cedrela_na_por_variavel <- cedrela_na_por_variavel[cedrela_na_por_variavel > 0]

# Ordenar da mais problemática para a menos problemática
sort(cedrela_na_por_variavel, decreasing = TRUE)


plot(
  stack_cedrela_final[[1]],
  main = "Pontos com e sem NA nas variáveis ambientais"
)

points(
  cedrela_pts$lon,
  cedrela_pts$lat,
  pch = 20,
  col = "blue",
  cex = 0.5
)

points(
  cedrela_pts_NA$lon,
  cedrela_pts_NA$lat,
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
cedrela_pts_modelo <- cedrela_pts[!cedrela_pontos_com_NA, ]

cat("Cedrela - pontos após thinning:", nrow(cedrela_pts), "\n")
cat("Cedrela - pontos usados no MaxEnt:", nrow(cedrela_pts_modelo), "\n")
cat("Cedrela - pontos removidos por NA:", sum(cedrela_pontos_com_NA), "\n")

# -------------------------------------------------------------------------
# 3.4 Ajuste MaxEnt para Cedrela
# -------------------------------------------------------------------------

# Configuração final recomendada:
# - validação cruzada com 10 folds
# - 1000 iterações máximas
# - saída cloglog
# - jackknife
# - curvas de resposta
# - regra Maximum training sensitivity plus specificity

cedrela_args_maxent <- c(
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

cat("\n[6/8] Ajustando MaxEnt para Cedrela...\n")
cat("Início:", format(Sys.time(), "%d/%m/%Y %H:%M:%S"), "\n")
tempo_cedrela <- system.time({

set.seed(1350)
cedrela_modelo <- dismo::maxent(
  x = stack_cedrela_final,
  p = cedrela_pts_modelo[, c("lon", "lat")],
  args = cedrela_args_maxent,
  path = out_cedrela,
  silent = FALSE
)

})

cat("Fim:", format(Sys.time(), "%d/%m/%Y %H:%M:%S"), "\n")
cat("Tempo total do ajuste MaxEnt - Cedrela:\n")
print(tempo_cedrela)

print(cedrela_modelo)

#--------------------------------------------------------------------------------------
# Gráfico de importancia de variáveis
# -------------------------------------------------------------------------------------

# Extrair resultados do modelo
cedrela_resultados <- cedrela_modelo@results

# Linhas com percentual de contribuição
linhas_contrib <- grep("\\.contribution$", rownames(cedrela_resultados), value = TRUE)

# Nomes das variáveis
vars <- gsub("\\.contribution$", "", linhas_contrib)

# Média da contribuição entre réplicas
contrib_media <- rowMeans(
  as.matrix(cedrela_resultados[linhas_contrib, , drop = FALSE]),
  na.rm = TRUE
)

# Linhas com importância por permutação
linhas_perm <- paste0(vars, ".permutation.importance")

perm_media <- rowMeans(
  as.matrix(cedrela_resultados[linhas_perm, , drop = FALSE]),
  na.rm = TRUE
)

# Montar tabela
imp_cedrela <- data.frame(
  variavel = vars,
  contribuicao = contrib_media,
  importancia_permutacao = perm_media
)

# Ordenar
imp_cedrela <- imp_cedrela[order(imp_cedrela$contribuicao, decreasing = TRUE), ]

imp_cedrela




ggplot(imp_cedrela,
       aes(x = reorder(variavel, contribuicao),
           y = contribuicao)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Variável ambiental",
    y = "Contribuição (%)",
    title = "Importância das variáveis - Cedrela",
    subtitle = "Percentual de contribuição no modelo MaxEnt"
  ) +
  theme_bw()



library(tidyr)
library(ggplot2)

imp_cedrela_long <- imp_cedrela |>
  pivot_longer(
    cols = c(contribuicao, importancia_permutacao),
    names_to = "metrica",
    values_to = "valor"
  )

ggplot(imp_cedrela_long,
       aes(x = reorder(variavel, valor),
           y = valor,
           fill = metrica)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    x = "Variável ambiental",
    y = "Importância (%)",
    fill = "Métrica",
    title = "Importância das variáveis - Cedrela"
  ) +
  theme_bw()


# -------------------------------------------------------------------------
# 3.5 AUC com validação cruzada para Cedrela
# -------------------------------------------------------------------------

# Verificar estrutura
rownames(cedrela_resultados)
colnames(cedrela_resultados)

# -------------------------------------------------------------------------
# Identificar colunas das réplicas da validação cruzada
# -------------------------------------------------------------------------

# Em geral, as colunas aparecem como:
# species_0, species_1, ..., species_9
# ou Cedrela_0, Cedrela_1, ..., Cedrela_9

cedrela_cols <- colnames(cedrela_resultados)

cedrela_cols_rep <- cedrela_cols[
  grepl("_[0-9]+$", cedrela_cols)
]

# Caso o padrão acima não funcione, remover colunas de média/desvio
if(length(cedrela_cols_rep) == 0){
  cedrela_cols_rep <- cedrela_cols[
    !grepl("average|avg|std|stddev|sd", cedrela_cols, ignore.case = TRUE)
  ]
}

cat("Colunas das réplicas usadas no cálculo da AUC:\n")
print(cedrela_cols_rep)

# -------------------------------------------------------------------------
# Identificar linhas exatas de AUC
# -------------------------------------------------------------------------

linha_auc_treino <- which(rownames(cedrela_resultados) == "Training.AUC")
linha_auc_teste  <- which(rownames(cedrela_resultados) == "Test.AUC")

if(length(linha_auc_treino) == 0){
  stop("Linha 'Training.AUC' não encontrada em cedrela_modelo@results.")
}

if(length(linha_auc_teste) == 0){
  stop("Linha 'Test.AUC' não encontrada em cedrela_modelo@results.")
}

# -------------------------------------------------------------------------
# Extrair AUC por fold
# -------------------------------------------------------------------------

cedrela_auc_cv <- data.frame(
  fold = seq_along(cedrela_cols_rep),
  replica = cedrela_cols_rep,
  AUC_treino = as.numeric(cedrela_resultados[linha_auc_treino, cedrela_cols_rep]),
  AUC_teste  = as.numeric(cedrela_resultados[linha_auc_teste,  cedrela_cols_rep])
)

# Diferença treino - teste como diagnóstico de possível sobreajuste
cedrela_auc_cv$delta_treino_teste <- cedrela_auc_cv$AUC_treino - cedrela_auc_cv$AUC_teste

print(cedrela_auc_cv)

# -------------------------------------------------------------------------
# Resumo da AUC por validação cruzada
# -------------------------------------------------------------------------

cedrela_auc_resumo <- data.frame(
  metrica = c("AUC treino", "AUC teste - validação cruzada", "Diferença treino - teste"),
  media = c(
    mean(cedrela_auc_cv$AUC_treino, na.rm = TRUE),
    mean(cedrela_auc_cv$AUC_teste, na.rm = TRUE),
    mean(cedrela_auc_cv$delta_treino_teste, na.rm = TRUE)
  ),
  desvio_padrao = c(
    sd(cedrela_auc_cv$AUC_treino, na.rm = TRUE),
    sd(cedrela_auc_cv$AUC_teste, na.rm = TRUE),
    sd(cedrela_auc_cv$delta_treino_teste, na.rm = TRUE)
  ),
  minimo = c(
    min(cedrela_auc_cv$AUC_treino, na.rm = TRUE),
    min(cedrela_auc_cv$AUC_teste, na.rm = TRUE),
    min(cedrela_auc_cv$delta_treino_teste, na.rm = TRUE)
  ),
  maximo = c(
    max(cedrela_auc_cv$AUC_treino, na.rm = TRUE),
    max(cedrela_auc_cv$AUC_teste, na.rm = TRUE),
    max(cedrela_auc_cv$delta_treino_teste, na.rm = TRUE)
  ),
  n_folds = c(
    sum(!is.na(cedrela_auc_cv$AUC_treino)),
    sum(!is.na(cedrela_auc_cv$AUC_teste)),
    sum(!is.na(cedrela_auc_cv$delta_treino_teste))
  )
)

print(cedrela_auc_resumo)

# -------------------------------------------------------------------------
# AUC principal a reportar no artigo
# -------------------------------------------------------------------------

cedrela_auc_treino_media <- mean(cedrela_auc_cv$AUC_treino, na.rm = TRUE)
cedrela_auc_treino_sd    <- sd(cedrela_auc_cv$AUC_treino, na.rm = TRUE)

cedrela_auc_teste_media <- mean(cedrela_auc_cv$AUC_teste, na.rm = TRUE)
cedrela_auc_teste_sd    <- sd(cedrela_auc_cv$AUC_teste, na.rm = TRUE)

cedrela_delta_auc_media <- mean(cedrela_auc_cv$delta_treino_teste, na.rm = TRUE)

cat("\nCedrela - AUC treino:",
    round(cedrela_auc_treino_media, 3),
    "±",
    round(cedrela_auc_treino_sd, 3),
    "\n")

cat("Cedrela - AUC teste por validação cruzada:",
    round(cedrela_auc_teste_media, 3),
    "±",
    round(cedrela_auc_teste_sd, 3),
    "\n")

cat("Cedrela - diferença média AUC treino - teste:",
    round(cedrela_delta_auc_media, 3),
    "\n")

# -------------------------------------------------------------------------
# 3.6 Predição espacial contínua para Cedrela
# -------------------------------------------------------------------------

cat("\nGerando predição espacial para Cedrela...\n")

cedrela_pred_bruto <- predict(
  cedrela_modelo,
  stack_cedrela_final,
  progress = "text"
)

class(cedrela_pred_bruto)

# Se o MaxEnt retornar uma camada por fold,
# calcular o mapa médio e o desvio-padrão entre folds
if(inherits(cedrela_pred_bruto, "RasterStack") | 
   inherits(cedrela_pred_bruto, "RasterBrick")){
  
  cedrela_pred <- raster::calc(
    cedrela_pred_bruto,
    fun = mean,
    na.rm = TRUE
  )
  
  cedrela_pred_sd <- raster::calc(
    cedrela_pred_bruto,
    fun = sd,
    na.rm = TRUE
  )
  
} else {
  
  cedrela_pred <- cedrela_pred_bruto
  cedrela_pred_sd <- NULL
  
}

# -------------------------------------------------------------------------
# Mapa contínuo médio em ggplot2
# -------------------------------------------------------------------------

library(ggplot2)
library(sf)

# Converter raster para data.frame
cedrela_pred_df <- raster::as.data.frame(
  cedrela_pred,
  xy = TRUE,
  na.rm = TRUE
)

# Padronizar nome da coluna de predição
names(cedrela_pred_df)[3] <- "adequabilidade"

# Converter contorno do Brasil para sf, se disponível
brasil_sf <- NULL

if(exists("brasil_uf") && !is.null(brasil_uf)){
  brasil_sf <- sf::st_as_sf(brasil_uf)
}

if(is.null(brasil_sf) && exists("brasil_sp") && !is.null(brasil_sp)){
  brasil_sf <- sf::st_as_sf(brasil_sp)
}

# Extensão do raster
ext_cedrela <- raster::extent(cedrela_pred)

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
mapa_cedrela_gg <- ggplot() +
  geom_raster(
    data = cedrela_pred_df,
    aes(
      x = x,
      y = y,
      fill = adequabilidade
    )
  ) +
  camada_brasil +
  geom_point(
    data = cedrela_pts_modelo,
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
    xlim = c(ext_cedrela@xmin, ext_cedrela@xmax),
    ylim = c(ext_cedrela@ymin, ext_cedrela@ymax),
    expand = FALSE
  ) +
  labs(
    title = expression(italic("Cedrela") ~ "- adequabilidade MaxEnt"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid = element_line(linewidth = 0.2, color = "gray85")
  )

mapa_cedrela_gg

#Salvar figura
ggsave(
  filename =  "figuras/cedrela_adequabilidade.png",
  plot = mapa_cedrela_gg,
  width = 8,
  height = 6,
  dpi = 300
)


# -------------------------------------------------------------------------
# Mapa de desvio-padrão entre folds em ggplot2
# -------------------------------------------------------------------------

if(!is.null(cedrela_pred_sd)){
  
  # Converter raster para data.frame
  cedrela_pred_sd_df <- raster::as.data.frame(
    cedrela_pred_sd,
    xy = TRUE,
    na.rm = TRUE
  )
  
  # Padronizar nome da coluna
  names(cedrela_pred_sd_df)[3] <- "desvio_padrao"
  
  # Converter contorno do Brasil para sf, se disponível
  brasil_sf <- NULL
  
  if(exists("brasil_uf") && !is.null(brasil_uf)){
    brasil_sf <- sf::st_as_sf(brasil_uf)
  }
  
  if(is.null(brasil_sf) && exists("brasil_sp") && !is.null(brasil_sp)){
    brasil_sf <- sf::st_as_sf(brasil_sp)
  }
  
  # Extensão do raster
  ext_cedrela_sd <- raster::extent(cedrela_pred_sd)
  
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
  mapa_cedrela_sd_gg <- ggplot() +
    geom_raster(
      data = cedrela_pred_sd_df,
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
      xlim = c(ext_cedrela_sd@xmin, ext_cedrela_sd@xmax),
      ylim = c(ext_cedrela_sd@ymin, ext_cedrela_sd@ymax),
      expand = FALSE
    ) +
    labs(
      title = expression(italic("Cedrela") ~ "- desvio-padrão entre folds"),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      panel.grid = element_line(linewidth = 0.2, color = "gray85")
    )
  
  mapa_cedrela_sd_gg
}


if(!is.null(cedrela_pred_sd)){
  
  ggsave(
    filename = "figuras/cedrela_desvio_padrao.png",
    plot = mapa_cedrela_sd_gg,
    width = 8,
    height = 6,
    dpi = 300
  )
  
}

# -------------------------------------------------------------------------
# 3.7 Limiar maxSSS e TSS para Cedrela
# -------------------------------------------------------------------------

# Verificar métricas de threshold disponíveis
grep(
  "Maximum.training.sensitivity.plus.specificity",
  rownames(cedrela_resultados),
  value = TRUE
)

# Buscar preferencialmente o limiar em Cloglog
cedrela_linhas_limiar <- grep(
  "Maximum.training.sensitivity.plus.specificity.Cloglog.threshold",
  rownames(cedrela_resultados)
)

# Caso o nome exato não exista, buscar qualquer threshold maxSSS
if(length(cedrela_linhas_limiar) == 0){
  cedrela_linhas_limiar <- grep(
    "Maximum.training.sensitivity.plus.specificity.*threshold",
    rownames(cedrela_resultados)
  )
}

if(length(cedrela_linhas_limiar) == 0){
  stop("Limiar Maximum training sensitivity plus specificity não encontrado.")
}

# Extrair limiares apenas das réplicas/folds
cedrela_limiar_vals <- as.numeric(
  cedrela_resultados[cedrela_linhas_limiar, cedrela_cols_rep, drop = FALSE]
)

cedrela_limiar_vals <- cedrela_limiar_vals[!is.na(cedrela_limiar_vals)]

cedrela_limiar_maxSSS <- mean(cedrela_limiar_vals, na.rm = TRUE)
cedrela_limiar_maxSSS_sd <- sd(cedrela_limiar_vals, na.rm = TRUE)

cat("Cedrela - limiar maxSSS médio:",
    round(cedrela_limiar_maxSSS, 3),
    "±",
    round(cedrela_limiar_maxSSS_sd, 3),
    "\n")

# -------------------------------------------------------------------------
# Avaliar sensibilidade, especificidade, TSS e omissão no limiar maxSSS
# -------------------------------------------------------------------------

set.seed(1350)

cedrela_bg_tss <- dismo::randomPoints(
  mask = stack_cedrela_final,
  n = 10000,
  ext = ext_brasil
)

colnames(cedrela_bg_tss) <- c("lon", "lat")
cedrela_bg_tss <- as.data.frame(cedrela_bg_tss)

# Extrair adequabilidade nos pontos de presença usados no modelo
cedrela_s_pres <- raster::extract(
  cedrela_pred,
  cedrela_pts_modelo[, c("lon", "lat")]
)

# Extrair adequabilidade nos pontos de background
cedrela_s_bg <- raster::extract(
  cedrela_pred,
  cedrela_bg_tss[, c("lon", "lat")]
)

# Remover NAs
cedrela_s_pres <- cedrela_s_pres[!is.na(cedrela_s_pres)]
cedrela_s_bg <- cedrela_s_bg[!is.na(cedrela_s_bg)]

# Matriz de confusão usando presença e background
cedrela_TP <- sum(cedrela_s_pres >= cedrela_limiar_maxSSS)
cedrela_FN <- sum(cedrela_s_pres <  cedrela_limiar_maxSSS)

cedrela_FP <- sum(cedrela_s_bg >= cedrela_limiar_maxSSS)
cedrela_TN <- sum(cedrela_s_bg <  cedrela_limiar_maxSSS)

# Métricas
cedrela_sensibilidade <- cedrela_TP / (cedrela_TP + cedrela_FN)
cedrela_especificidade <- cedrela_TN / (cedrela_TN + cedrela_FP)

cedrela_TSS <- cedrela_sensibilidade + cedrela_especificidade - 1
cedrela_omissao <- cedrela_FN / (cedrela_TP + cedrela_FN)

cedrela_metricas_threshold <- data.frame(
  especie = "Cedrela",
  limiar_maxSSS = cedrela_limiar_maxSSS,
  sensibilidade = cedrela_sensibilidade,
  especificidade = cedrela_especificidade,
  TSS = cedrela_TSS,
  omissao = cedrela_omissao,
  TP = cedrela_TP,
  FN = cedrela_FN,
  FP = cedrela_FP,
  TN = cedrela_TN
)

print(cedrela_metricas_threshold)

cat("Cedrela - sensibilidade:", round(cedrela_sensibilidade, 3), "\n")
cat("Cedrela - especificidade:", round(cedrela_especificidade, 3), "\n")
cat("Cedrela - TSS:", round(cedrela_TSS, 3), "\n")
cat("Cedrela - omissão:", round(cedrela_omissao, 3), "\n")


# -------------------------------------------------------------------------
# 3.8 Mapas classificados para Cedrela em ggplot2
# -------------------------------------------------------------------------

library(ggplot2)
library(sf)
library(dplyr)

# -------------------------------------------------------------------------
# Reclassificação: 4 classes de adequabilidade
# -------------------------------------------------------------------------

cedrela_rcl_4classes <- matrix(
  c(
    -Inf, 0.265, 0,
    0.265, 0.530, 1,
    0.530, 0.765, 2,
    0.765, Inf,   3
  ),
  ncol = 3,
  byrow = TRUE
)

cedrela_mapa4 <- raster::reclassify(
  cedrela_pred,
  rcl = cedrela_rcl_4classes
)

# -------------------------------------------------------------------------
# Reclassificação: mapa binário maxSSS
# -------------------------------------------------------------------------

cedrela_rcl_bin <- matrix(
  c(
    -Inf, cedrela_limiar_maxSSS, 0,
    cedrela_limiar_maxSSS, Inf, 1
  ),
  ncol = 3,
  byrow = TRUE
)

cedrela_mapa_bin <- raster::reclassify(
  cedrela_pred,
  rcl = cedrela_rcl_bin
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
ext_cedrela <- raster::extent(cedrela_pred)

# -------------------------------------------------------------------------
# Mapa 4 classes em ggplot2
# -------------------------------------------------------------------------

cedrela_mapa4_df <- raster::as.data.frame(
  cedrela_mapa4,
  xy = TRUE,
  na.rm = TRUE
)

names(cedrela_mapa4_df)[3] <- "classe_valor"

cedrela_mapa4_df <- cedrela_mapa4_df |>
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

mapa_cedrela_4classes_gg <- ggplot() +
  geom_raster(
    data = cedrela_mapa4_df,
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
    xlim = c(ext_cedrela@xmin, ext_cedrela@xmax),
    ylim = c(ext_cedrela@ymin, ext_cedrela@ymax),
    expand = FALSE
  ) +
  labs(
    title = expression(italic("Cedrela") ~ "- classes de adequabilidade"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid = element_line(linewidth = 0.2, color = "gray85")
  )

mapa_cedrela_4classes_gg

#Salvar
ggsave(
  filename = "figuras/cedrela_4_classes_ggplot.png",
  plot = mapa_cedrela_4classes_gg,
  width = 8,
  height = 6,
  dpi = 300
)


# -------------------------------------------------------------------------
# Mapa binário em ggplot2
# -------------------------------------------------------------------------

cedrela_mapa_bin_df <- raster::as.data.frame(
  cedrela_mapa_bin,
  xy = TRUE,
  na.rm = TRUE
)

names(cedrela_mapa_bin_df)[3] <- "classe_valor"

cedrela_mapa_bin_df <- cedrela_mapa_bin_df |>
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

mapa_cedrela_binario_gg <- ggplot() +
  geom_raster(
    data = cedrela_mapa_bin_df,
    aes(
      x = x,
      y = y,
      fill = classe
    )
  ) +
  camada_brasil +
  geom_point(
    data = cedrela_pts_modelo,
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
    xlim = c(ext_cedrela@xmin, ext_cedrela@xmax),
    ylim = c(ext_cedrela@ymin, ext_cedrela@ymax),
    expand = FALSE
  ) +
  labs(
    title = paste0(
      "Cedrela - binário maxSSS = ",
      round(cedrela_limiar_maxSSS, 3)
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

mapa_cedrela_binario_gg



#Salvar figura
ggsave(
  filename = "figuras/cedrela_binario_maxSSS_ggplot.png",
  plot = mapa_cedrela_binario_gg,
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------------
# 3.9 Área por classe para Cedrela
# -------------------------------------------------------------------------

# Área de cada célula em km²
cedrela_area_pixel_4 <- raster::area(cedrela_mapa4, na.rm = TRUE)

cedrela_vals_classe4 <- raster::values(cedrela_mapa4)
cedrela_vals_area4 <- raster::values(cedrela_area_pixel_4)

cedrela_ok4 <- !is.na(cedrela_vals_classe4) & !is.na(cedrela_vals_area4)

cedrela_vals_classe4 <- cedrela_vals_classe4[cedrela_ok4]
cedrela_vals_area4 <- cedrela_vals_area4[cedrela_ok4]

cedrela_area4_soma <- tapply(
  cedrela_vals_area4,
  cedrela_vals_classe4,
  sum
)

cedrela_area4_total <- sum(cedrela_area4_soma)

cedrela_areas4 <- data.frame(
  valor = as.integer(names(cedrela_area4_soma)),
  classe = c("Inadequado", "Baixa", "Média", "Alta")[
    as.integer(names(cedrela_area4_soma)) + 1
  ],
  area_km2 = round(as.numeric(cedrela_area4_soma), 0),
  pct = round(
    100 * as.numeric(cedrela_area4_soma) / cedrela_area4_total,
    2
  )
)

print(cedrela_areas4)

# Área binária
cedrela_area_pixel_bin <- raster::area(cedrela_mapa_bin, na.rm = TRUE)

cedrela_vals_bin <- raster::values(cedrela_mapa_bin)
cedrela_vals_area_bin <- raster::values(cedrela_area_pixel_bin)

cedrela_ok_bin <- !is.na(cedrela_vals_bin) & !is.na(cedrela_vals_area_bin)

cedrela_vals_bin <- cedrela_vals_bin[cedrela_ok_bin]
cedrela_vals_area_bin <- cedrela_vals_area_bin[cedrela_ok_bin]

cedrela_area_bin_soma <- tapply(
  cedrela_vals_area_bin,
  cedrela_vals_bin,
  sum
)

cedrela_area_bin_total <- sum(cedrela_area_bin_soma)

cedrela_areasbin <- data.frame(
  valor = as.integer(names(cedrela_area_bin_soma)),
  classe = c("Inadequado", "Adequado")[
    as.integer(names(cedrela_area_bin_soma)) + 1
  ],
  area_km2 = round(as.numeric(cedrela_area_bin_soma), 0),
  pct = round(
    100 * as.numeric(cedrela_area_bin_soma) / cedrela_area_bin_total,
    2
  )
)

print(cedrela_areasbin)


# -------------------------------------------------------------------------
# 3.10 Exportação dos resultados de Cedrela
# -------------------------------------------------------------------------

# Raster contínuo médio
raster::writeRaster(
  cedrela_pred,
  file.path(out_dir, "cedrela_adequabilidade_continua_media_cv.tif"),
  overwrite = TRUE
)

# Raster de desvio-padrão entre folds
if(!is.null(cedrela_pred_sd)){
  
  raster::writeRaster(
    cedrela_pred_sd,
    file.path(out_dir, "cedrela_adequabilidade_sd_cv.tif"),
    overwrite = TRUE
  )
  
}

# Raster com 4 classes
raster::writeRaster(
  cedrela_mapa4,
  file.path(out_dir, "cedrela_4_classes.tif"),
  overwrite = TRUE,
  datatype = "INT1U"
)

# Raster binário maxSSS
raster::writeRaster(
  cedrela_mapa_bin,
  file.path(out_dir, "cedrela_binario_maxSSS.tif"),
  overwrite = TRUE,
  datatype = "INT1U"
)

# Tabelas de resultados
writexl::write_xlsx(
  list(
    AUC_por_fold = cedrela_auc_cv,
    AUC_resumo = cedrela_auc_resumo,
    Metricas_threshold = cedrela_metricas_threshold,
    Area_4_classes = cedrela_areas4,
    Area_binaria = cedrela_areasbin
  ),
  path = file.path(out_dir, "cedrela_resultados_maxent_cv.xlsx")
)

# Objetos R
saveRDS(
  cedrela_modelo,
  file.path(out_dir, "modelo_cedrela_cv.rds")
)

saveRDS(
  cedrela_pred,
  file.path(out_dir, "predicao_cedrela_media_cv.rds")
)

if(!is.null(cedrela_pred_sd)){
  saveRDS(
    cedrela_pred_sd,
    file.path(out_dir, "predicao_cedrela_sd_cv.rds")
  )
}

# Salvar ambiente de trabalho
save.image("data/modelo_cedrela_cv.RData")
