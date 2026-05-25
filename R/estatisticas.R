# Instalar pacotes (apenas na primeira vez)
install.packages("sf")
install.packages("dplyr")

# Carregar bibliotecas
library(sf)
library(dplyr)

# =========================
# LER CSVs COMBINADOS
# =========================
# Gerados pelo script plot_pontos.R

handroanthus_combinado <- read.csv("Handroanthus_combinado.csv")

cedrela_combinado <- read.csv("Cedrela_combinado.csv")

# =========================
# CONVERTER PARA OBJETO ESPACIAL
# =========================
# CRS 4326 (graus). st_distance usa cálculo geodésico
# automaticamente e retorna distâncias em metros.

handroanthus_sf <- st_as_sf(
  handroanthus_combinado,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

cedrela_sf <- st_as_sf(
  cedrela_combinado,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

# =========================
# FUNÇÃO DE ESTATÍSTICAS
# =========================
# Calcula a matriz de distâncias par-a-par e retorna
# min, média, mediana, sd e máxima. As distâncias são
# convertidas para km. A diagonal (distância de cada
# ponto consigo mesmo = 0) é removida.

estatisticas_distancia <- function(pontos_sf, rotulo) {
  
  # Matriz de distâncias par-a-par (em metros)
  d <- st_distance(pontos_sf)
  
  # Converter para numérico e passar para km
  d <- as.numeric(d) / 1000
  
  # Reorganizar como matriz
  n <- nrow(pontos_sf)
  d <- matrix(d, nrow = n, ncol = n)
  
  # Descartar diagonal (zeros) usando apenas o
  # triângulo superior, sem a diagonal
  d_pares <- d[upper.tri(d, diag = FALSE)]
  
  # Estatísticas
  cat("=========================================\n")
  cat(rotulo, "(", n, "pontos)\n")
  cat("Total de pares:", length(d_pares), "\n")
  cat("Distância mínima (km):", min(d_pares), "\n")
  cat("Distância média   (km):", mean(d_pares), "\n")
  cat("Distância mediana (km):", median(d_pares), "\n")
  cat("Desvio padrão     (km):", sd(d_pares), "\n")
  cat("Distância máxima  (km):", max(d_pares), "\n")
  cat("=========================================\n\n")
  
  # Retornar como data.frame para uso posterior
  data.frame(
    dataset = rotulo,
    n_pontos = n,
    n_pares = length(d_pares),
    min_km = min(d_pares),
    media_km = mean(d_pares),
    mediana_km = median(d_pares),
    sd_km = sd(d_pares),
    max_km = max(d_pares)
  )
}

# =========================
# HANDROANTHUS
# =========================

stats_handroanthus <- estatisticas_distancia(
  handroanthus_sf,
  "Handroanthus combinado"
)

# =========================
# CEDRELA
# =========================

stats_cedrela <- estatisticas_distancia(
  cedrela_sf,
  "Cedrela combinado"
)

# =========================
# SALVAR RESULTADOS
# =========================

stats_final <- bind_rows(
  stats_handroanthus,
  stats_cedrela
)

write.csv(
  stats_final,
  "estatisticas_distancia.csv",
  row.names = FALSE
)

cat("Arquivo estatisticas_distancia.csv salvo.\n\n")

# =========================
# PREVISÃO DE THINNING
# =========================
# Para cada distância mínima d (em km), simula um
# thinning guloso: percorre os pontos e mantém apenas
# aqueles que estão a >= d km de TODOS os pontos já
# mantidos. Retorna quantos pontos sobrariam.
#
# Observação: a ordem dos pontos afeta o resultado
# exato, mas a contagem final tende a ser estável.
# Usamos a ordem original do CSV.

prever_thinning <- function(pontos_sf, distancias_km, rotulo) {
  
  n <- nrow(pontos_sf)
  
  # Matriz de distâncias par-a-par (em km)
  d <- st_distance(pontos_sf)
  d <- matrix(as.numeric(d) / 1000, nrow = n, ncol = n)
  
  resultados <- data.frame()
  
  cat("=========================================\n")
  cat("Previsão de thinning -", rotulo, "\n")
  cat("Total inicial:", n, "pontos\n")
  cat("-----------------------------------------\n")
  
  for (dist_min in distancias_km) {
    
    # Vetor lógico: TRUE = ponto mantido
    manter <- rep(TRUE, n)
    
    for (i in seq_len(n)) {
      if (!manter[i]) next
      # Marca como descartado qualquer ponto j > i
      # que esteja a menos de dist_min de i
      vizinhos <- which(d[i, ] < dist_min & manter)
      vizinhos <- vizinhos[vizinhos != i]
      manter[vizinhos] <- FALSE
    }
    
    n_mantidos <- sum(manter)
    n_removidos <- n - n_mantidos
    pct_mantidos <- round(100 * n_mantidos / n, 1)
    
    cat(
      "Thinning ", dist_min, " km: ",
      n_mantidos, " mantidos | ",
      n_removidos, " removidos | ",
      pct_mantidos, "%\n",
      sep = ""
    )
    
    resultados <- bind_rows(
      resultados,
      data.frame(
        dataset = rotulo,
        thinning_km = dist_min,
        n_inicial = n,
        n_mantidos = n_mantidos,
        n_removidos = n_removidos,
        pct_mantidos = pct_mantidos
      )
    )
  }
  
  cat("=========================================\n\n")
  
  resultados
}

# Distâncias a testar
distancias_thinning <- c(1, 5, 10, 15, 20)

# Handroanthus
thin_handroanthus <- prever_thinning(
  handroanthus_sf,
  distancias_thinning,
  "Handroanthus combinado"
)

# Cedrela
thin_cedrela <- prever_thinning(
  cedrela_sf,
  distancias_thinning,
  "Cedrela combinado"
)

# Salvar resultados de thinning
thin_final <- bind_rows(
  thin_handroanthus,
  thin_cedrela
)

write.csv(
  thin_final,
  "previsao_thinning.csv",
  row.names = FALSE
)

cat("Arquivo previsao_thinning.csv salvo.\n\n")

# =========================
# RELATÓRIO EM TXT
# =========================
# Gera um arquivo de texto com descrição das métricas
# e os resultados das estatísticas de distância e
# da previsão de thinning para cada dataset.

# Função auxiliar: formata uma linha de stats
formatar_stats <- function(s) {
  paste0(
    "  Número de pontos:    ", s$n_pontos, "\n",
    "  Número de pares:     ", s$n_pares, "\n",
    "  Distância mínima:    ", round(s$min_km, 3), " km\n",
    "  Distância média:     ", round(s$media_km, 3), " km\n",
    "  Distância mediana:   ", round(s$mediana_km, 3), " km\n",
    "  Desvio padrão (sd):  ", round(s$sd_km, 3), " km\n",
    "  Distância máxima:    ", round(s$max_km, 3), " km\n"
  )
}

# Função auxiliar: formata o bloco de thinning
formatar_thinning <- function(thin_df) {
  linhas <- ""
  for (i in seq_len(nrow(thin_df))) {
    r <- thin_df[i, ]
    linhas <- paste0(
      linhas,
      "  Thinning ", sprintf("%2d", r$thinning_km), " km: ",
      r$n_mantidos, " pontos mantidos de ", r$n_inicial,
      " (", r$pct_mantidos, "%) | ",
      r$n_removidos, " removidos\n"
    )
  }
  linhas
}

relatorio <- paste0(
  "===========================================================\n",
  "  RELATÓRIO DE ESTATÍSTICAS DE DISTÂNCIA E THINNING\n",
  "  Gerado em: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "===========================================================\n\n",
  
  "DATASETS ANALISADOS\n",
  "-----------------------------------------------------------\n",
  "Os dados combinam pontos de ocorrência da planilha de\n",
  "isótopos (filtrados pela Amazônia Legal) com pontos de\n",
  "ocorrência do Brasil com variáveis ambientais, removendo\n",
  "duplicatas por lon/lat. Foram analisados dois gêneros:\n",
  "Handroanthus e Cedrela.\n\n",
  
  "MÉTRICAS DE DISTÂNCIA\n",
  "-----------------------------------------------------------\n",
  "As distâncias foram calculadas par-a-par entre todos os\n",
  "pontos de cada dataset, usando a função st_distance() do\n",
  "pacote sf com CRS geográfico (EPSG:4326). Isso garante\n",
  "cálculo geodésico (distância sobre a superfície da Terra),\n",
  "com resultado em metros, convertido para km. Apenas pares\n",
  "únicos (triângulo superior da matriz, sem a diagonal)\n",
  "foram considerados.\n\n",
  "  - Distância mínima: menor distância entre quaisquer dois\n",
  "    pontos. Valores muito baixos indicam pontos quase\n",
  "    coincidentes, sinal de possível redundância amostral.\n",
  "  - Distância média: média aritmética de todas as\n",
  "    distâncias par-a-par. Sensível a outliers (pontos\n",
  "    muito distantes do restante).\n",
  "  - Distância mediana: valor central das distâncias,\n",
  "    robusto a outliers. Reflete melhor a distância típica.\n",
  "  - Desvio padrão (sd): dispersão das distâncias em torno\n",
  "    da média. Valor alto indica grande variabilidade\n",
  "    espacial entre as ocorrências.\n",
  "  - Distância máxima: maior separação entre dois pontos\n",
  "    do dataset, refletindo a extensão geográfica total.\n\n",
  
  "PREVISÃO DE THINNING\n",
  "-----------------------------------------------------------\n",
  "O thinning espacial é uma técnica usada em modelagem de\n",
  "distribuição de espécies (SDM) para reduzir o viés\n",
  "amostral causado por concentração de pontos em regiões\n",
  "mais amostradas. Consiste em manter apenas pontos\n",
  "separados por uma distância mínima d.\n\n",
  "O algoritmo aplicado é guloso: percorre os pontos na\n",
  "ordem do arquivo e, para cada ponto ainda mantido,\n",
  "descarta todos os demais que estiverem a menos de d km\n",
  "dele. O número final de pontos é razoavelmente estável,\n",
  "mas os pontos exatos que sobram dependem da ordem.\n\n",
  "Foram testadas cinco distâncias: 1, 5, 10, 15 e 20 km.\n\n",
  
  "===========================================================\n",
  "  RESULTADOS - HANDROANTHUS COMBINADO\n",
  "===========================================================\n\n",
  "Estatísticas de distância (km):\n",
  formatar_stats(stats_handroanthus), "\n",
  "Previsão de thinning:\n",
  formatar_thinning(thin_handroanthus), "\n",
  
  "===========================================================\n",
  "  RESULTADOS - CEDRELA COMBINADO\n",
  "===========================================================\n\n",
  "Estatísticas de distância (km):\n",
  formatar_stats(stats_cedrela), "\n",
  "Previsão de thinning:\n",
  formatar_thinning(thin_cedrela), "\n",
  
  "===========================================================\n",
  "  ARQUIVOS GERADOS\n",
  "===========================================================\n",
  "  - estatisticas_distancia.csv  (métricas em formato CSV)\n",
  "  - previsao_thinning.csv       (previsão de thinning em CSV)\n",
  "  - relatorio_distancias.txt    (este relatório)\n"
)

writeLines(
  relatorio,
  "relatorio_distancias.txt"
)

cat("Arquivo relatorio_distancias.txt salvo.\n")