# Instalar pacotes (apenas na primeira vez)
install.packages("sf")
install.packages("ggplot2")
install.packages("dplyr")

# Carregar bibliotecas
library(sf)
library(ggplot2)
library(dplyr)

# =========================
# LER SHAPEFILE AMAZÔNIA
# =========================

amazonia <- st_read("shapefile/amazonia_legal.shp")

# =========================
# LER CSVs
# =========================

handroanthus <- read.csv("Handroanthus_filtrado.csv")

cedrela <- read.csv("Cedrela_filtrado.csv")

handroanthus_br <- read.csv("Handroanthus_BR_filtrado.csv")

cedrela_br <- read.csv("Cedrela_BR_filtrado.csv")

# =========================
# CONVERTER PARA OBJETO ESPACIAL
# =========================

handroanthus_sf <- st_as_sf(
  handroanthus,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

cedrela_sf <- st_as_sf(
  cedrela,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

handroanthus_br_sf <- st_as_sf(
  handroanthus_br,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

cedrela_br_sf <- st_as_sf(
  cedrela_br,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

# =========================
# MAPA HANDROANTHUS
# =========================

ggplot() +
  
  # Contorno Amazônia Legal
  geom_sf(
    data = amazonia,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Pontos
  geom_sf(
    data = handroanthus_sf,
    color = "blue",
    size = 1.5,
    alpha = 0.7
  ) +
  
  ggtitle("Handroanthus - Amazônia Legal") +
  
  theme_bw()

# Salvar figura
ggsave(
  "mapa_handroanthus.png",
  width = 8,
  height = 8,
  dpi = 300
)

# =========================
# MAPA CEDRELA
# =========================

ggplot() +
  
  # Contorno Amazônia Legal
  geom_sf(
    data = amazonia,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Pontos
  geom_sf(
    data = cedrela_sf,
    color = "red",
    size = 1.5,
    alpha = 0.7
  ) +
  
  ggtitle("Cedrela - Amazônia Legal") +
  
  theme_bw()

# Salvar figura
ggsave(
  "mapa_cedrela.png",
  width = 8,
  height = 8,
  dpi = 300
)

# =========================
# MAPA HANDROANTHUS BR
# =========================

ggplot() +
  
  # Contorno Amazônia Legal
  geom_sf(
    data = amazonia,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Pontos
  geom_sf(
    data = handroanthus_br_sf,
    color = "darkgreen",
    size = 1.2,
    alpha = 0.6
  ) +
  
  ggtitle("Handroanthus BR - Amazônia Legal") +
  
  theme_bw()

# Salvar figura
ggsave(
  "mapa_handroanthus_br.png",
  width = 8,
  height = 8,
  dpi = 300
)

# =========================
# MAPA CEDRELA BR
# =========================

ggplot() +
  
  # Contorno Amazônia Legal
  geom_sf(
    data = amazonia,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Pontos
  geom_sf(
    data = cedrela_br_sf,
    color = "darkorange",
    size = 1.2,
    alpha = 0.6
  ) +
  
  ggtitle("Cedrela BR - Amazônia Legal") +
  
  theme_bw()

# Salvar figura
ggsave(
  "mapa_cedrela_br.png",
  width = 8,
  height = 8,
  dpi = 300
)

# =========================
# COMBINAR DATASETS POR GÊNERO
# =========================
# Junta os dois conjuntos de Handroanthus (planilha + BR var_amb)
# e os dois de Cedrela, mantendo apenas as colunas comuns
# (lat, long, x, y) para permitir o bind.

handroanthus_combinado <- bind_rows(
  handroanthus %>% select(lat, long, x, y),
  handroanthus_br %>% select(lat, long, x, y)
)

cedrela_combinado <- bind_rows(
  cedrela %>% select(lat, long, x, y),
  cedrela_br %>% select(lat, long, x, y)
)

# Remover duplicatas por lon/lat após a junção
handroanthus_combinado <- handroanthus_combinado %>%
  distinct(long, lat, .keep_all = TRUE)

cedrela_combinado <- cedrela_combinado %>%
  distinct(long, lat, .keep_all = TRUE)

# Salvar CSVs combinados
write.csv(
  handroanthus_combinado,
  "Handroanthus_combinado.csv",
  row.names = FALSE
)

write.csv(
  cedrela_combinado,
  "Cedrela_combinado.csv",
  row.names = FALSE
)

# Converter para objeto espacial
handroanthus_combinado_sf <- st_as_sf(
  handroanthus_combinado,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

cedrela_combinado_sf <- st_as_sf(
  cedrela_combinado,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

# Contar pontos
cat(
  "Handroanthus combinado:",
  nrow(handroanthus_combinado_sf),
  "pontos\n"
)

cat(
  "Cedrela combinado:",
  nrow(cedrela_combinado_sf),
  "pontos\n"
)

# =========================
# MAPA HANDROANTHUS COMBINADO
# =========================

ggplot() +
  
  # Contorno Amazônia Legal
  geom_sf(
    data = amazonia,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Pontos
  geom_sf(
    data = handroanthus_combinado_sf,
    color = "blue",
    size = 1.2,
    alpha = 0.6
  ) +
  
  ggtitle("Handroanthus combinado - Amazônia Legal") +
  
  theme_bw()

# Salvar figura
ggsave(
  "mapa_handroanthus_combinado.png",
  width = 8,
  height = 8,
  dpi = 300
)

# =========================
# MAPA CEDRELA COMBINADO
# =========================

ggplot() +
  
  # Contorno Amazônia Legal
  geom_sf(
    data = amazonia,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Pontos
  geom_sf(
    data = cedrela_combinado_sf,
    color = "red",
    size = 1.2,
    alpha = 0.6
  ) +
  
  ggtitle("Cedrela combinado - Amazônia Legal") +
  
  theme_bw()

# Salvar figura
ggsave(
  "mapa_cedrela_combinado.png",
  width = 8,
  height = 8,
  dpi = 300
)