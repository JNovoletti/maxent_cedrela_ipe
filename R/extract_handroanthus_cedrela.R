# Instalar pacotes (apenas na primeira vez)
install.packages("readxl")
install.packages("dplyr")
install.packages("sf")

# Carregar bibliotecas
library(readxl)
library(dplyr)
library(sf)

# Ler arquivo Excel
dados <- read_excel("C:/Users/Aluno/Downloads/+Planilha_madeiras_14_04_26_compilada.xlsx", sheet="Resultado_isotopos")


# =========================
# LER SHAPEFILE
# =========================
# Coloque todos os arquivos do shapefile
# (.shp, .dbf, .shx etc.) na mesma pasta

amazonia <- st_read("shapefile/amazonia_legal.shp")

# =========================
# HANDROANTHUS
# =========================

handroanthus <- dados %>%
  filter(Genus == "Handroanthus")

# Remover duplicatas por latitude/longitude
handroanthus_unico <- handroanthus %>%
  distinct(latitude, longitude, .keep_all = TRUE)

# Criar colunas lat, long, x e y
handroanthus_unico <- handroanthus_unico %>%
  mutate(
    lat = latitude,
    long = longitude,
    x = longitude,
    y = latitude
  )

# Selecionar colunas desejadas
handroanthus_final <- handroanthus_unico %>%
  select(
    lat,
    long,
    x,
    y,
    Species,
    Scientific_name,
    Genus,
    Code,
    Cod.Lab,
    Family,
    Point,
    Site,
    State,
    UC
  )

# Converter para objeto espacial
handroanthus_sf <- st_as_sf(
  handroanthus_final,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

# Filtrar apenas pontos dentro da Amazônia Legal
handroanthus_amazonia <- st_intersection(
  handroanthus_sf,
  amazonia
)

# Contar pontos
cat(
  "Handroanthus:",
  nrow(handroanthus_amazonia),
  "pontos na Amazônia Legal\n"
)

# Remover geometria
handroanthus_saida <- st_drop_geometry(
  handroanthus_amazonia
)

# Salvar CSV
write.csv(
  handroanthus_saida,
  "Handroanthus_filtrado.csv",
  row.names = FALSE
)

# =========================
# CEDRELA
# =========================

cedrela <- dados %>%
  filter(Genus == "Cedrela")

# Remover duplicatas por latitude/longitude
cedrela_unico <- cedrela %>%
  distinct(latitude, longitude, .keep_all = TRUE)

# Criar colunas lat, long, x e y
cedrela_unico <- cedrela_unico %>%
  mutate(
    lat = latitude,
    long = longitude,
    x = longitude,
    y = latitude
  )

# Selecionar colunas desejadas
cedrela_final <- cedrela_unico %>%
  select(
    lat,
    long,
    x,
    y,
    Species,
    Scientific_name,
    Genus,
    Code,
    Cod.Lab,
    Family,
    Point,
    Site,
    State,
    UC
  )

# Converter para objeto espacial
cedrela_sf <- st_as_sf(
  cedrela_final,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

# Filtrar apenas pontos dentro da Amazônia Legal
cedrela_amazonia <- st_intersection(
  cedrela_sf,
  amazonia
)

# Contar pontos
cat(
  "Cedrela:",
  nrow(cedrela_amazonia),
  "pontos na Amazônia Legal\n"
)

# Remover geometria
cedrela_saida <- st_drop_geometry(
  cedrela_amazonia
)

# Salvar CSV
write.csv(
  cedrela_saida,
  "Cedrela_filtrado.csv",
  row.names = FALSE
)

# =========================
# HANDROANTHUS BR (var_amb)
# =========================
# Lê arquivo de pontos do Brasil inteiro com variáveis ambientais.
# As colunas ambientais NÃO são preservadas; mantemos apenas lon/lat/x/y.

handroanthus_br <- read.csv(
  "data/handroanthus_var_amb.csv"
)

# Manter apenas colunas de coordenadas
handroanthus_br <- handroanthus_br %>%
  select(lon, lat, x, y)

# Remover duplicatas por lon/lat
handroanthus_br_unico <- handroanthus_br %>%
  distinct(lon, lat, .keep_all = TRUE)

# Padronizar nome (long) para casar com o restante do script
handroanthus_br_final <- handroanthus_br_unico %>%
  mutate(long = lon) %>%
  select(lat, long, x, y)

# Converter para objeto espacial
handroanthus_br_sf <- st_as_sf(
  handroanthus_br_final,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

# Filtrar apenas pontos dentro da Amazônia Legal
handroanthus_br_amazonia <- st_intersection(
  handroanthus_br_sf,
  amazonia
)

# Contar pontos
cat(
  "Handroanthus BR:",
  nrow(handroanthus_br_amazonia),
  "pontos na Amazônia Legal\n"
)

# Remover geometria
handroanthus_br_saida <- st_drop_geometry(
  handroanthus_br_amazonia
)

# Salvar CSV
write.csv(
  handroanthus_br_saida,
  "Handroanthus_BR_filtrado.csv",
  row.names = FALSE
)

# =========================
# CEDRELA BR (var_amb)
# =========================
# Lê arquivo de pontos do Brasil inteiro com variáveis ambientais.
# As colunas ambientais NÃO são preservadas; mantemos apenas lon/lat/x/y.

cedrela_br <- read.csv(
  "data/cedrela_br_var_amb.csv"
)

# Manter apenas colunas de coordenadas
cedrela_br <- cedrela_br %>%
  select(lon, lat, x, y)

# Remover duplicatas por lon/lat
cedrela_br_unico <- cedrela_br %>%
  distinct(lon, lat, .keep_all = TRUE)

# Padronizar nome (long) para casar com o restante do script
cedrela_br_final <- cedrela_br_unico %>%
  mutate(long = lon) %>%
  select(lat, long, x, y)

# Converter para objeto espacial
cedrela_br_sf <- st_as_sf(
  cedrela_br_final,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

# Filtrar apenas pontos dentro da Amazônia Legal
cedrela_br_amazonia <- st_intersection(
  cedrela_br_sf,
  amazonia
)

# Contar pontos
cat(
  "Cedrela BR:",
  nrow(cedrela_br_amazonia),
  "pontos na Amazônia Legal\n"
)

# Remover geometria
cedrela_br_saida <- st_drop_geometry(
  cedrela_br_amazonia
)

# Salvar CSV
write.csv(
  cedrela_br_saida,
  "Cedrela_BR_filtrado.csv",
  row.names = FALSE
)

# =========================
# FINAL
# =========================

cat("Arquivos CSV gerados com sucesso!\n")