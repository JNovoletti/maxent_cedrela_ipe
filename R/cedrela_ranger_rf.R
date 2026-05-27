# =========================================================================
# SDM Cedrela – Random Forest (ranger) + VSURF
# Área de estudo  : Amazônia Legal
# Fontes de dados : GBIF Brasil (data/cedrela_br_var_amb.csv)
#                   Amostras de campo (Planilha_madeiras_14_04_26_compilada.xlsx)
# Rasters         : rasters/1.*.tif … 76.*.tif  → máscara Amazônia
# Thinning        : 5 km (spThin)
# =========================================================================

# rm(list = ls()); gc(reset = TRUE); graphics.off()

# =========================================================================
# 1. PACOTES
# =========================================================================

pkgs <- c(
  "readxl", "writexl", "tidyverse", "terra", "sf",
  "here", "geobr", "spThin", "ranger", "VSURF", "ggplot2", "pROC"
)
invisible(lapply(pkgs, function(p) {
  if (!require(p, character.only = TRUE))
    install.packages(p, dependencies = TRUE)
  library(p, character.only = TRUE)
}))

set.seed(1350)
RNGkind("Mersenne-Twister", "Inversion", "Rounding")

out_dir <- here::here("output_ranger", "cedrela")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(here::here("figuras"))) dir.create(here::here("figuras"))


# =========================================================================
# 2. SHAPEFILE DA AMAZÔNIA LEGAL
# =========================================================================

shp_path <- here::here("shapefile", "amazonia_legal.shp")

if (file.exists(shp_path)) {
  amazonia <- sf::st_read(shp_path, quiet = TRUE)
} else {
  cat("Shapefile local não encontrado. Baixando via geobr::read_amazon()...\n")
  amazonia <- geobr::read_amazon(year = 2012, showProgress = FALSE)
  dir.create(here::here("shapefile"), showWarnings = FALSE)
  sf::st_write(amazonia, shp_path, quiet = TRUE)
  cat("Shapefile salvo em:", shp_path, "\n")
}

amazonia <- sf::st_make_valid(amazonia)
if (is.na(sf::st_crs(amazonia)) || sf::st_crs(amazonia)$epsg != 4326)
  amazonia <- sf::st_transform(amazonia, 4326)


# =========================================================================
# 3. RASTERS → STACK E MÁSCARA PARA AMAZÔNIA
# =========================================================================

raster_files <- list.files(
  here::here("rasters"),
  pattern = "\\.tif$",
  full.names = TRUE
)

excluir <- "Brazil_env_stack|Amazon_stack|Caatinga_elev_corridor|elev_corridor_minus1|48\\.density|17\\.fossilfuels_dep"
raster_files <- raster_files[!grepl(excluir, basename(raster_files))]

r_brasil <- terra::rast(raster_files)

nomes <- tools::file_path_sans_ext(basename(raster_files))
nomes <- gsub("^[0-9]+\\.", "", nomes)
names(r_brasil) <- nomes

amz_vect <- terra::vect(amazonia)
amz_vect <- terra::project(amz_vect, terra::crs(r_brasil))

r_amazonia <- terra::crop(r_brasil, amz_vect)
r_amazonia <- terra::mask(r_amazonia, amz_vect)

cat("Camadas ambientais:", terra::nlyr(r_amazonia), "\n")
cat("Nomes:", paste(names(r_amazonia), collapse = ", "), "\n")


# =========================================================================
# 4. OCORRÊNCIAS
# =========================================================================

# --- 4.1 GBIF / Brasil (CSV com variáveis ambientais) --------------------

gbif_br <- tryCatch(
  read.csv(here::here("data", "cedrela_br_var_amb.csv")),
  error = function(e) {
    message("cedrela_br_var_amb.csv não encontrado; usando cedrela_br.csv")
    read.csv(here::here("data", "cedrela_br.csv"))
  }
)

lon_col <- intersect(c("lon", "longitude", "x", "decimalLongitude"), names(gbif_br))[1]
lat_col <- intersect(c("lat", "latitude",  "y", "decimalLatitude"),  names(gbif_br))[1]

gbif_pts <- data.frame(
  lon   = as.numeric(gbif_br[[lon_col]]),
  lat   = as.numeric(gbif_br[[lat_col]]),
  fonte = "gbif_br"
)

# --- 4.2 Amostras de campo (Excel) ---------------------------------------

excel_path <- "C:/Users/Aluno/Downloads/+Planilha_madeiras_14_04_26_compilada.xlsx"

amostras_raw <- tryCatch(
  readxl::read_excel(excel_path, sheet = "Resultado_isotopos"),
  error = function(e) {
    message("Excel não encontrado ou aba diferente: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(amostras_raw)) {
  amostras_cedrela <- amostras_raw %>%
    dplyr::filter(Genus == "Cedrela") %>%
    dplyr::mutate(
      lon = as.numeric(longitude),
      lat = as.numeric(latitude)
    ) %>%
    dplyr::select(lon, lat) %>%
    dplyr::mutate(fonte = "campo")
} else {
  amostras_cedrela <- data.frame(lon = numeric(0), lat = numeric(0), fonte = character(0))
}

# --- 4.3 Combinar, remover NA e duplicatas -------------------------------

occ_all <- dplyr::bind_rows(gbif_pts, amostras_cedrela) %>%
  dplyr::filter(
    !is.na(lon), !is.na(lat),
    lon >= -180, lon <= 180,
    lat >= -90,  lat <= 90
  ) %>%
  dplyr::distinct(lon, lat, .keep_all = TRUE)

cat("Ocorrências totais (antes do clip):", nrow(occ_all), "\n")

# --- 4.4 Clip para Amazônia Legal ----------------------------------------

occ_sf  <- sf::st_as_sf(occ_all, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
occ_amz <- sf::st_filter(occ_sf, amazonia)
occ_pts <- sf::st_drop_geometry(occ_amz)

cat("Ocorrências dentro da Amazônia:", nrow(occ_pts), "\n")


# =========================================================================
# 5. THINNING ESPACIAL DE 5 KM
# =========================================================================

occ_thin_in <- occ_pts %>%
  dplyr::select(lon, lat) %>%
  dplyr::mutate(especie = "Cedrela")

set.seed(1350)
thin_list <- spThin::thin(
  loc.data               = occ_thin_in,
  lat.col                = "lat",
  long.col               = "lon",
  spec.col               = "especie",
  thin.par               = 5,
  reps                   = 50,
  locs.thinned.list.return = TRUE,
  write.files            = FALSE,
  write.log.file         = FALSE,
  verbose                = FALSE
)

melhor_rep <- which.max(sapply(thin_list, nrow))
occ_final  <- thin_list[[melhor_rep]]
names(occ_final) <- c("lon", "lat")

cat("Cedrela – pontos após thinning de 5 km:", nrow(occ_final), "\n")


# =========================================================================
# 6. PSEUDO-AUSÊNCIAS (BACKGROUND) DENTRO DA AMAZÔNIA
# =========================================================================

set.seed(1350)
bg_pts <- terra::spatSample(
  r_amazonia[[1]],
  size   = 1000,
  method = "random",
  na.rm  = TRUE,
  as.df  = FALSE,
  xy     = TRUE
)

bg_df <- as.data.frame(bg_pts)[, c("x", "y")]
names(bg_df) <- c("lon", "lat")

cat("Pontos de background gerados:", nrow(bg_df), "\n")


# =========================================================================
# 7. EXTRAÇÃO DE VARIÁVEIS AMBIENTAIS
# =========================================================================

pres_sv  <- terra::vect(occ_final, geom = c("lon", "lat"), crs = "EPSG:4326")
pres_sv  <- terra::project(pres_sv, terra::crs(r_amazonia))
pres_env <- terra::extract(r_amazonia, pres_sv, method = "bilinear", ID = FALSE)
pres_env$presence <- 1L

bg_sv   <- terra::vect(bg_df, geom = c("lon", "lat"), crs = "EPSG:4326")
bg_sv   <- terra::project(bg_sv, terra::crs(r_amazonia))
bg_env  <- terra::extract(r_amazonia, bg_sv, method = "bilinear", ID = FALSE)
bg_env$presence <- 0L

model_data <- dplyr::bind_rows(pres_env, bg_env)
model_data  <- model_data[complete.cases(model_data), ]

cat("Linhas totais para modelagem:", nrow(model_data),
    "| presença:", sum(model_data$presence),
    "| background:", sum(model_data$presence == 0), "\n")

resposta   <- as.factor(model_data$presence)
preditoras <- model_data %>% dplyr::select(-presence)

if (anyNA(preditoras)) {
  cat("Aviso: valores ausentes em preditoras – removidos.\n")
  ok <- complete.cases(preditoras)
  preditoras <- preditoras[ok, ]
  resposta   <- resposta[ok]
}


# =========================================================================
# 8. SELEÇÃO DE VARIÁVEIS – VSURF
# =========================================================================

cat("\nRodando VSURF para seleção de variáveis...\n")
set.seed(1350)

vsurf_res <- VSURF::VSURF(
  x          = as.matrix(preditoras),
  y          = resposta,
  ntree      = 500,
  nfor.thres = 20,
  nfor.interp= 100,
  nfor.pred  = 10,
  nsd        = 1,
  parallel   = FALSE,
  verbose    = FALSE
)

vars_pred <- names(preditoras)[vsurf_res$varselect.interp]

cat("Variáveis selecionadas (varselect.interp):\n"); print(vars_pred)
cat("Total selecionado:", length(vars_pred), "\n")

preditoras_sel <- preditoras[, vars_pred, drop = FALSE]
model_data_sel <- cbind(preditoras_sel, presence = resposta)


# =========================================================================
# 9. MODELO RANDOM FOREST – ranger (probability = TRUE)
# =========================================================================

formula_rf <- as.formula(paste("presence ~", paste(vars_pred, collapse = " + ")))

cat("\nAjustando Random Forest (ranger)...\n")
set.seed(1350)

rf_final <- ranger::ranger(
  formula      = formula_rf,
  data         = model_data_sel,
  num.trees    = 500,
  probability  = TRUE,
  importance   = "impurity",
  min.node.size= 10,
  seed         = 1350
)

cat("OOB prediction error:", round(rf_final$prediction.error, 4), "\n")

imp_df <- data.frame(
  variavel    = names(rf_final$variable.importance),
  importancia = as.numeric(rf_final$variable.importance)
) %>%
  dplyr::arrange(dplyr::desc(importancia))

print(imp_df)

gg_imp <- ggplot(imp_df, aes(x = reorder(variavel, importancia), y = importancia)) +
  geom_col(fill = "#d73027") +
  coord_flip() +
  labs(
    title = expression(italic("Cedrela") ~ "– Importância das variáveis (RF)"),
    x     = "Variável",
    y     = "Importância (impurity)"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(here::here("figuras"), "cedrela_rf_importancia.png"),
  plot = gg_imp, width = 8, height = 6, dpi = 300
)


# =========================================================================
# 10. MÉTRICAS DE AVALIAÇÃO (AUC, maxSSS, TSS)
# =========================================================================

# OOB predictions (ranger probability = TRUE armazena em $predictions)
oob_probs  <- rf_final$predictions[, "1"]
oob_labels <- as.integer(as.character(model_data_sel$presence))

# AUC via curva ROC
roc_obj <- pROC::roc(oob_labels, oob_probs, quiet = TRUE)
auc_val <- as.numeric(pROC::auc(roc_obj))

# maxSSS: limiar que maximiza Sensibilidade + Especificidade
coords_roc <- pROC::coords(
  roc_obj, x = "all", input = "threshold",
  ret       = c("threshold", "sensitivity", "specificity"),
  transpose = FALSE
)
coords_roc$tss <- coords_roc$sensitivity + coords_roc$specificity - 1
idx_max        <- which.max(coords_roc$tss)

limiar_maxSSS  <- coords_roc$threshold[idx_max]
sensibilidade  <- coords_roc$sensitivity[idx_max]
especificidade <- coords_roc$specificity[idx_max]
TSS            <- coords_roc$tss[idx_max]
omissao        <- 1 - sensibilidade

# Matriz de confusão no limiar maxSSS
pred_bin <- as.integer(oob_probs >= limiar_maxSSS)
TP <- sum(pred_bin == 1L & oob_labels == 1L)
FN <- sum(pred_bin == 0L & oob_labels == 1L)
FP <- sum(pred_bin == 1L & oob_labels == 0L)
TN <- sum(pred_bin == 0L & oob_labels == 0L)

acuracia <- (TP + TN) / (TP + TN + FP + FN)
precisao <- TP / (TP + FP)
f1       <- 2 * TP / (2 * TP + FP + FN)
n_tot    <- TP + TN + FP + FN
pe       <- ((TP + FP) / n_tot) * ((TP + FN) / n_tot) +
            ((FN + TN) / n_tot) * ((FP + TN) / n_tot)
kappa    <- (acuracia - pe) / (1 - pe)

metricas_rf <- data.frame(
  metrica = c(
    "AUC_OOB", "OOB_erro_predicao",
    "Limiar_maxSSS", "TSS", "Sensibilidade", "Especificidade",
    "Omissao", "Acuracia", "Precisao", "F1_score", "Kappa",
    "TP", "FN", "FP", "TN"
  ),
  valor = c(
    round(auc_val, 4), round(rf_final$prediction.error, 4),
    round(limiar_maxSSS, 4), round(TSS, 4), round(sensibilidade, 4), round(especificidade, 4),
    round(omissao, 4), round(acuracia, 4), round(precisao, 4), round(f1, 4), round(kappa, 4),
    TP, FN, FP, TN
  )
)

cat("\n--- Métricas de avaliação (OOB) ---\n")
print(metricas_rf)

cat("\nAUC (OOB):    ", round(auc_val, 4), "\n")
cat("Limiar maxSSS:", round(limiar_maxSSS, 4), "\n")
cat("TSS:          ", round(TSS, 4), "\n")
cat("Sensibilidade:", round(sensibilidade, 4), "\n")
cat("Especificidade:", round(especificidade, 4), "\n")
cat("Omissão:      ", round(omissao, 4), "\n")


# =========================================================================
# 11. MAPA DE PROBABILIDADE DE OCORRÊNCIA
# =========================================================================

r_sel <- r_amazonia[[vars_pred]]

cat("\nGerando mapa de probabilidade de ocorrência...\n")

pred_fun <- function(model, data, ...) {
  pred <- predict(model, data = as.data.frame(data))
  as.numeric(pred$predictions[, "1"])
}

prob_map <- terra::predict(r_sel, rf_final, fun = pred_fun, na.rm = TRUE)
names(prob_map) <- "prob_ocorrencia"

# --- Mapa contínuo -------------------------------------------------------

prob_df <- as.data.frame(prob_map, xy = TRUE, na.rm = TRUE)
names(prob_df)[3] <- "probabilidade"

brasil_uf <- tryCatch(
  geobr::read_state(year = 2020, showProgress = FALSE),
  error = function(e) NULL
)
brasil_sf <- if (!is.null(brasil_uf)) sf::st_as_sf(brasil_uf) else NULL

mapa_gg <- ggplot() +
  geom_raster(data = prob_df, aes(x = x, y = y, fill = probabilidade)) +
  { if (!is.null(brasil_sf))
      geom_sf(data = brasil_sf, fill = NA, color = "gray30",
              linewidth = 0.3, inherit.aes = FALSE) } +
  geom_point(
    data = occ_final, aes(x = lon, y = lat),
    shape = 21, fill = "yellow", color = "black", size = 1, stroke = 0.2
  ) +
  scale_fill_viridis_c(
    name     = "Probabilidade\nde ocorrência",
    limits   = c(0, 1),
    na.value = "transparent"
  ) +
  labs(
    title = expression(italic("Cedrela") ~ "– probabilidade de ocorrência (RF)"),
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid      = element_line(linewidth = 0.2, color = "gray85")
  )

mapa_gg

ggsave(
  file.path(here::here("figuras"), "cedrela_probabilidade_ocorrencia.png"),
  plot = mapa_gg, width = 8, height = 6, dpi = 300
)

# --- Mapa binário (limiar maxSSS) ----------------------------------------

bin_map <- terra::classify(
  prob_map,
  matrix(c(-Inf, limiar_maxSSS, 0, limiar_maxSSS, Inf, 1), ncol = 3, byrow = TRUE)
)
names(bin_map) <- "presenca_binario"

bin_df <- as.data.frame(bin_map, xy = TRUE, na.rm = TRUE) %>%
  dplyr::mutate(classe = factor(presenca_binario, levels = c(0, 1),
                                labels = c("Ausente", "Presente")))

mapa_bin_gg <- ggplot() +
  geom_raster(data = bin_df, aes(x = x, y = y, fill = classe)) +
  { if (!is.null(brasil_sf))
      geom_sf(data = brasil_sf, fill = NA, color = "gray30",
              linewidth = 0.3, inherit.aes = FALSE) } +
  geom_point(
    data = occ_final, aes(x = lon, y = lat),
    shape = 21, fill = "yellow", color = "black", size = 1, stroke = 0.2
  ) +
  scale_fill_manual(
    name   = "Classe",
    values = c("Ausente" = "#d9d9d9", "Presente" = "#a50026"),
    drop   = FALSE, na.value = "transparent"
  ) +
  labs(
    title = paste0("Cedrela – presença binária (maxSSS = ", round(limiar_maxSSS, 3), ")"),
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid      = element_line(linewidth = 0.2, color = "gray85")
  )

mapa_bin_gg

ggsave(
  file.path(here::here("figuras"), "cedrela_binario_rf.png"),
  plot = mapa_bin_gg, width = 8, height = 6, dpi = 300
)


# =========================================================================
# 12. EXPORTAÇÃO
# =========================================================================

terra::writeRaster(
  prob_map,
  file.path(out_dir, "cedrela_prob_ocorrencia.tif"),
  overwrite = TRUE
)

terra::writeRaster(
  bin_map,
  file.path(out_dir, "cedrela_binario_rf.tif"),
  overwrite = TRUE, datatype = "INT1U"
)

writexl::write_xlsx(
  list(
    variaveis_selecionadas = data.frame(variavel = vars_pred),
    importancia_variaveis  = imp_df,
    metricas_avaliacao     = metricas_rf,
    ocorrencias_thinning   = occ_final
  ),
  path = file.path(out_dir, "cedrela_rf_resultados.xlsx")
)

saveRDS(rf_final,  file.path(out_dir, "modelo_cedrela_ranger.rds"))
saveRDS(vsurf_res, file.path(out_dir, "vsurf_cedrela.rds"))

save.image(file.path(out_dir, "cedrela_ranger_workspace.RData"))

cat("\n=== Cedrela RF concluído ===\n")
cat("Outputs salvos em:", out_dir, "\n")
