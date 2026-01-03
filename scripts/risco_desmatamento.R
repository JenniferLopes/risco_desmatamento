# ANÁLISE DE DESMATAMENTO E PRODUÇÃO AGRÍCOLA NO BRASIL (2015-2021)

# Projeto: Análise da relação entre produção agrícola e desmatamento
# Objetivo: Identificar padrões espaciais e temporais entre produção agrícola
#           (soja, milho, arroz) e desmatamento usando dados do IBGE PAM e INPE PRODES

# Fontes de dados:
# - PAM (Produção Agrícola Municipal) - IBGE via Base dos Dados
# - PRODES (Monitoramento de Desmatamento) - INPE via Base dos Dados

# 1. Configuração inicial ----------------------------------------------------

# install.packages("pacman")
library(pacman)

pacman::p_load(
  dplyr,
  tidyr,
  ggplot2,
  scales,
  zoo,
  DBI,
  dbplyr,
  bigrquery,
  sf,
  geobr,
  leaflet,
  htmlwidgets,
  randomForest,
  skimr,
  rsample,
  yardstick,
  tibble,
  stringr,
  readr,
  htmltools,
  usethis)

# Paleta de cores para mapas
CORES_MAPA <- list(
  verde_floresta = "#2E7D32",
  verde_seco    = "#8D9B6A",
  amarelo_queimado = "#E6B800",
  marrom_terra  = "#8D4B2D",
  vermelho_escuro = "#7F1D1D")

# Paleta de cores para gráficos
CORES_GRAFICO <- list(
  azul_medio  = "#3A5FA8",
  azul_claro  = "#8FAADC",
  marrom_claro = "#A07878",
  bege_medio  = "#E6DED6",
  cinza_quente = "#6E6E6E")

# Tema para ggplot2
tema_custom <- function() {
  ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0),
      plot.subtitle = ggplot2::element_text(
        size = 12,
        color = CORES_GRAFICO$cinza_quente),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "#E5E5E5",
        linewidth = 0.3),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", size = 11),
      axis.title = ggplot2::element_text(face = "bold", size = 12),
      plot.margin = ggplot2::margin(
        t = 20, r = 20, b = 20, l = 20,
        unit = "pt"))
}

# 2. Conectar ao BigQuery e extrair dados PAM e PRODES -----------------------

# Conexão com BigQuery
con <- DBI::dbConnect(
  bigrquery::bigquery(),
  project = "desmatamento-482501",
  use_legacy_sql = FALSE,
  quiet = TRUE)

# # Validação da conexão
# teste_conexao <- DBI::dbGetQuery(con, "SELECT 1 AS test")
# stopifnot(teste_conexao$test == 1)

# Extração PAM
pam_tbl <- tbl(con, "basedosdados.br_ibge_pam.lavoura_temporaria")

dados_pam <- pam_tbl |>
  filter(
    ano >= 2015, ano <= 2021,
    produto %in% c("Soja (em grão)", "Milho (em grão)", "Arroz (em casca)"),
    quantidade_produzida > 0,
    area_colhida > 0) |>
  group_by(ano, id_municipio, sigla_uf) |>
  summarise(
    producao_total_ton = sum(quantidade_produzida, na.rm = TRUE),
    area_colhida_ha    = sum(area_colhida, na.rm = TRUE),
    valor_producao_mil = sum(valor_producao, na.rm = TRUE) / 1000,
    .groups = "drop") |>
  collect() |>
  mutate(
    id_municipio = as.character(id_municipio),
    ano         = as.integer(ano),
    produtividade = producao_total_ton / area_colhida_ha,
    produtividade = if_else(
      is.infinite(produtividade) | is.nan(produtividade),
      NA_real_,
      produtividade))

view(dados_pam)

# Extração PRODES
prodes_tbl <- tbl(con, "basedosdados.br_inpe_prodes.municipio_bioma")

dados_prodes <- prodes_tbl |>
  filter(ano >= 2015, ano <= 2021) |>
  group_by(ano, id_municipio) |>
  summarise(
    desmatamento_ha = sum(desmatado, na.rm = TRUE),
    area_total_km2  = sum(area_total, na.rm = TRUE),
    vegetacao_natural_km2 = sum(vegetacao_natural, na.rm = TRUE),
    .groups = "drop") |>
  collect() |>
  mutate(
    id_municipio = as.character(id_municipio),
    ano        = as.integer(ano),
    pct_vegetacao = (vegetacao_natural_km2 / area_total_km2) * 100,
    pct_vegetacao = if_else(
      is.infinite(pct_vegetacao) | is.nan(pct_vegetacao),
      0,
      pct_vegetacao))

view(dados_prodes)
# Fechar conexão
# DBI::dbDisconnect(con)

# 3. Integração e engenharia de features ----------------------------------

# Integração das bases PAM e PRODES
df_integrado <- dados_pam |>
  inner_join(dados_prodes, by = c("ano", "id_municipio")) |>
  arrange(id_municipio, ano)

view(df_integrado)

# Criação de features temporais e derivadas
df_features_temp <- df_integrado |>
  group_by(id_municipio) |>
  mutate(
    desmat_lag1 = lag(desmatamento_ha),
    desmat_lag2 = lag(desmatamento_ha, 2),
    desmat_ma3  = zoo::rollmeanr(desmatamento_ha, 3, fill = NA),
    crescimento_producao =
      (producao_total_ton - lag(producao_total_ton)) /
      (lag(producao_total_ton) + 1),
    expansao_area = area_colhida_ha - lag(area_colhida_ha),
    risco_desmat_ton =
      desmatamento_ha / (producao_total_ton + 1),
    pressao_vegetacao =
      area_colhida_ha / (vegetacao_natural_km2 * 100 + 1)) |>
  ungroup()

view(df_features_temp)

# Limpeza e consistência
df_limpo <- df_features_temp |>
  filter(
    !is.na(desmat_lag1),
    desmatamento_ha >= 0,
    producao_total_ton > 0,
    is.finite(risco_desmat_ton))

view(df_limpo)

# Definição das classes de risco
q <- quantile(df_limpo$risco_desmat_ton, c(0.33, 0.67), na.rm = TRUE)

df_completo <- df_limpo |>
  mutate(
    classe_risco = case_when(
      risco_desmat_ton < q[1] ~ "Baixo",
      risco_desmat_ton < q[2] ~ "Medio",
      TRUE ~ "Alto"),
    classe_risco = factor(classe_risco, levels = c("Baixo", "Medio", "Alto")))

view(df_completo)

# 4. Análise descritiva ---------------------------------------------------

# Visão geral utilizando skimr
skim_geral <- skimr::skim(df_completo)
view(skim_geral)

# Estatísticas por classe de risco
analise_por_risco <- df_completo |>
  group_by(classe_risco) |>
  summarise(
    n_municipios     = n_distinct(id_municipio),
    n_observacoes    = n(),
    desmat_medio     = mean(desmatamento_ha, na.rm = TRUE),
    desmat_mediano   = median(desmatamento_ha, na.rm = TRUE),
    producao_media   = mean(producao_total_ton, na.rm = TRUE),
    producao_mediana = median(producao_total_ton, na.rm = TRUE),
    area_media       = mean(area_colhida_ha, na.rm = TRUE),
    produtividade_media = mean(produtividade, na.rm = TRUE),
    .groups = "drop")

view(analise_por_risco)

# Estatísticas temporais agregadas
analise_temporal <- df_completo |>
  group_by(ano) |>
  summarise(
    desmatamento_total = sum(desmatamento_ha, na.rm = TRUE),
    producao_total     = sum(producao_total_ton, na.rm = TRUE),
    n_municipios       = n_distinct(id_municipio),
    .groups = "drop")

view(analise_temporal)

analise_descritiva <- list(
  visao_geral     = skim_geral,
  por_classe_risco = analise_por_risco,
  temporal        = analise_temporal)

# 5. Divisão de dados e modelagem --------------------------------------------

# Divisão temporal
dados_treino <- df_completo |> filter(ano < 2020)
dados_validacao <- df_completo |> filter(ano == 2020)
dados_teste <- df_completo |> filter(ano >= 2021)

splits <- list(
  treino    = dados_treino,
  validacao = dados_validacao,
  teste     = dados_teste)

# Definição das features para o modelo
FEATURES <- c(
  "producao_total_ton",
  "area_colhida_ha",
  "produtividade",
  "pct_vegetacao",
  "desmat_lag1",
  "desmat_lag2",
  "crescimento_producao",
  "expansao_area",
  "pressao_vegetacao")

# Preparação de dados de treino
dados_modelo_treino <- dados_treino |>
  select(desmatamento_ha, all_of(FEATURES)) |>
  na.omit()
view(dados_modelo_treino)

# Modelo baseline: regressão linear
modelo_lm <- lm(desmatamento_ha ~ ., data = dados_modelo_treino)

# Modelo Random Forest
modelo_rf <- randomForest(
  desmatamento_ha ~ .,
  data = dados_modelo_treino,
  ntree = 300,
  importance = TRUE,
  na.action = na.omit)

# Preparação de dados de teste para avaliação
dados_modelo_teste <- dados_teste |>
  select(id_municipio, ano, sigla_uf, classe_risco, desmatamento_ha, all_of(FEATURES)) |>
  na.omit()
view(dados_modelo_teste)

# Geração de previsões
predicoes <- predict(modelo_rf, newdata = dados_modelo_teste)

dados_preditos <- dados_modelo_teste |>
  mutate(pred = predicoes)

view(dados_preditos)

# Métricas padrão com yardstick
metricas_modelo <- dados_preditos |>
  yardstick::metrics(
    truth = desmatamento_ha,
    estimate = pred)

view(metricas_modelo)

# Métricas adicionais
metricas_complementares <- dados_preditos |>
  mutate(
    erro_abs = abs(pred - desmatamento_ha),
    erro_pct = (erro_abs / (desmatamento_ha + 1)) * 100 ) |>
  summarise(
    MAE = round(mean(erro_abs, na.rm = TRUE), 2),
    RMSE = round(
      sqrt(mean((pred - desmatamento_ha)^2, na.rm = TRUE)),
      2),
    R2 = round(
      cor(pred, desmatamento_ha, use = "complete.obs")^2,
      3),
    MAPE = round(mean(erro_pct, na.rm = TRUE), 2))

view(metricas_complementares)

avaliacao_modelo <- list(
  metricas_yardstick = metricas_modelo,
  metricas_resumo   = metricas_complementares,
  resultados        = dados_preditos)


# 6. Visualizações -------------------------------------------------------

# Criar diretórios de saída para gráficos
if (!dir.exists("outputs")) dir.create("outputs")
if (!dir.exists("outputs/figures")) dir.create("outputs/figures")


# Gráfico 1: Importância das variáveis ------------------------------------

imp_df <- importance(modelo_rf) |>
  as.data.frame() |>
  tibble::rownames_to_column("Feature") |>
  arrange(desc(`%IncMSE`)) |>
  slice_head(n = 10) |>
  mutate(
    Feature = stringr::str_replace_all(Feature, "_", " "),
    Feature = stringr::str_to_title(Feature))

view(imp_df)

grafico_importancia <- ggplot(imp_df,
                              aes(x = reorder(Feature, `%IncMSE`), y = `%IncMSE`)) +
  geom_col(fill = CORES_GRAFICO$azul_medio, alpha = 0.9, width = 0.7) +
  coord_flip() +
  labs(
    title    = "Importância das variáveis no modelo Random Forest",
    subtitle = "Variáveis por % de aumento no MSE.",
    x = NULL,
    y = "% de aumento no MSE", caption = "Jennifer Luz Lopes, 2025.") +
  tema_custom() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11))

grafico_importancia

ggsave("outputs/figures/grafico_1_importancia.png",
       grafico_importancia, width = 10, height = 6,
       dpi = 300, bg = "white")


# Gráfico 2: Scatter produção x desmatamento ------------------------------

df_filtrado <- df_completo |>
  filter(ano >= 2019)

n_amostra <- min(2000, nrow(df_filtrado))

df_plot_scatter <- df_filtrado |>
  slice_sample(n = n_amostra)

view(df_plot_scatter)

grafico_scatter <- ggplot(
  df_plot_scatter,
  aes(x = producao_total_ton / 1000, y = desmatamento_ha)) +
  geom_point(aes(color = classe_risco),
             alpha = 0.6, size = 2.5) +
  scale_color_manual(
    values = c(
      "Baixo" = CORES_MAPA$verde_floresta,
      "Medio" = CORES_MAPA$amarelo_queimado,
      "Alto"  = CORES_MAPA$vermelho_escuro),
    name = "Classe de risco") +
  scale_x_log10(
    labels = label_comma(),
    breaks = c(1, 10, 100, 1000)) +
  scale_y_log10(
    labels = label_comma(),
    breaks = c(10, 100, 1000, 10000)) +
  labs(
    title = "Relação entre produção agrícola e desmatamento",
    subtitle = sprintf(
      "Municípios brasileiros (2019-%d, n=%d) - cores representam risco relativo (ha/ton)",
      max(df_plot_scatter$ano), nrow(df_plot_scatter)),
    x = "Produção total (mil toneladas, escala log)",
    y = "Desmatamento (hectares, escala log)", caption = "Jennifer Luz Lopes, 2025.") +
  tema_custom() +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.9)))

grafico_scatter

ggsave("outputs/figures/grafico_2_scatter.png",
       grafico_scatter, width = 12, height = 7,
       dpi = 300, bg = "white")


# Gráfico 3: Série temporal -----------------------------------------------

# df_norm <- analise_descritiva$temporal |>
#   mutate(
#     desmat_z = scale(desmatamento_total)[,1],
#     prod_z   = scale(producao_total)[,1]
#   ) |>
#   select(ano, desmat_z, prod_z) |>
#   pivot_longer(
#     cols = c(desmat_z, prod_z),
#     names_to = "variavel",
#     values_to = "valor"
#   ) |>
#   mutate(
#     variavel = recode(
#       variavel,
#       desmat_z = "Desmatamento",
#       prod_z   = "Produção agrícola"
#     )
#   )
# 
# 
# 
# view(df_norm)
# 
# grafico_temporal <- ggplot(df_norm, aes(x = ano, y = valor, color = variavel)) +
#   geom_line(linewidth = 1.2) +
#   geom_point(size = 3) +
#   scale_color_manual(
#     values = c(
#       "Desmatamento"     = CORES_MAPA$vermelho_escuro,
#       "Produção agrícola" = CORES_GRAFICO$azul_medio),
#     name = NULL) +
#   scale_x_continuous(breaks = seq(2015, 2021, 1)) +
#   scale_y_continuous(
#     labels = label_percent(scale = 1),
#     limits = c(0, 105)) +
#   labs(
#     title    = "Evolução Temporal: Desmatamento vs Produção Agrícola",
#     subtitle = "Valores normalizados (0-100%)",
#     x        = "Ano",
#     y        = "Índice Normalizado") +
#   tema_custom() +
#   theme(legend.position = c(0.2, 0.9))
# 
# grafico_temporal
# 
# ggsave("outputs/figures/grafico_3_temporal.png",
#        grafico_temporal, width = 10, height = 6,
#        dpi = 300, bg = "white")


# 7. Mapas ----------------------------------------------------------------

# Criar diretórios de saída para mapas
if (!dir.exists("outputs")) dir.create("outputs")
if (!dir.exists("outputs/maps")) dir.create("outputs/maps")

# Definições
ano_ref <- 2021
n_top  <- 50

# Mapa 1: Municípios com maior desmatamento -------------------------------

top_municipios <- df_completo |>
  filter(ano == ano_ref) |>
  arrange(desc(desmatamento_ha)) |>
  slice_head(n = n_top)

view(top_municipios)

municipios_sf <- suppressMessages(
  geobr::read_municipality(year = 2020, showProgress = FALSE)) |>
  mutate(id_municipio = as.character(code_muni)) |>
  filter(id_municipio %in% top_municipios$id_municipio) |>
  left_join(top_municipios, by = "id_municipio")

view(municipios_sf)

pal_risco <- colorFactor(
  palette = c(
    CORES_MAPA$verde_floresta,
    CORES_MAPA$amarelo_queimado,
    CORES_MAPA$vermelho_escuro),
  domain = c("Baixo", "Medio", "Alto"))

mapa_top_risco <- leaflet(municipios_sf) |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(
    fillColor = ~pal_risco(classe_risco),
    fillOpacity = 0.8,
    color = "#FFFFFF",
    weight = 1.5,
    opacity = 1,
    label = ~lapply(paste0(
      "<strong style='font-size:14px;'>", name_muni, " - ", sigla_uf, "</strong><br/>",
      "<span style='font-size:12px;'>",
      "Desmatamento: <strong>", format(round(desmatamento_ha, 0), big.mark = "."), " ha</strong><br/>",
      "Produção: ", format(round(producao_total_ton / 1000, 1), big.mark = "."), " mil ton<br/>",
      "Classe: <strong style='color:",
      ifelse(classe_risco == "Alto", CORES_MAPA$vermelho_escuro,
             ifelse(classe_risco == "Medio", CORES_MAPA$amarelo_queimado,
                    CORES_MAPA$verde_floresta)),
      ";'>", classe_risco, "</strong>",
      "</span>"), htmltools::HTML),
    highlightOptions = highlightOptions(
      weight = 3,
      color = "#000000",
      fillOpacity = 0.9,
      bringToFront = TRUE)) |>
  addLegend(
    position = "bottomright",
    pal = pal_risco,
    values = ~classe_risco,
    title = "Classe de Risco",
    opacity = 1) |>
  addControl(
    html = paste0(
      "<div style='background:white; padding:15px; border-radius:8px; ",
      "box-shadow: 0 2px 10px rgba(0,0,0,0.2); font-family:Arial, sans-serif;'>",
      "<h3 style='margin:0 0 10px 0; color:#333; font-size:16px;'>",
     " Municípios de Maior Desmatamento</h3>",
      "<p style='margin:0; color:#666; font-size:13px;'>Ano: ", ano_ref, "</p>",
      "</div>"),
    position = "topleft")

htmlwidgets::saveWidget(
  mapa_top_risco,
  "outputs/maps/mapa_1_top_risco.html",
  selfcontained = TRUE)


# Mapa 1.1 estático

# Seleção dosmunicípios 
top_municipios <- df_completo |>
  filter(ano == ano_ref) |>
  arrange(desc(desmatamento_ha)) |>
  slice_head(n = n_top) |>
  mutate(id_municipio = as.character(id_municipio))

# Base espacial – TODOS os municípios do Brasil
municipios_sf <- suppressMessages(
  geobr::read_municipality(year = 2020, showProgress = FALSE)) |>
  mutate(id_municipio = as.character(code_muni))

# Join dos dados (SEM FILTRAR o mapa)
municipios_sf <- municipios_sf |>
  left_join(top_municipios, by = "id_municipio")

# Mapa estático
mapa_top_risco_static <- ggplot(municipios_sf) +
  geom_sf(
    aes(fill = classe_risco),
    color = "white",
    linewidth = 0.1,
    na.rm = FALSE) +
  scale_fill_manual(
    values = c(
      "Baixo" = CORES_MAPA$verde_floresta,
      "Medio" = CORES_MAPA$amarelo_queimado,
      "Alto"  = CORES_MAPA$vermelho_escuro),
    na.value = "grey90",
    name = "Classe de risco") +
  labs(
    title = "Municípios com maior desmatamento",
    subtitle = paste0("Top ", n_top, " municípios – ", ano_ref),
    caption = "Fonte: IBGE (PAM) e INPE (PRODES)") +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank())

# Plot
mapa_top_risco_static

# Mapa 2: Intensidade de desmatamento (coroplético) -----------------------

percentil_min <- 0.70

dados_mapa2 <- df_completo |>
  filter(ano == ano_ref) |>
  group_by(id_municipio, sigla_uf) |>
  summarise(
    desmatamento_ha = sum(desmatamento_ha, na.rm = TRUE),
    producao_total_ton = sum(producao_total_ton, na.rm = TRUE),
    risco_desmat_ton = mean(risco_desmat_ton, na.rm = TRUE),
    .groups = "drop") |>
  filter(desmatamento_ha > quantile(desmatamento_ha, percentil_min))

municipios_sf2 <- suppressMessages(
  geobr::read_municipality(year = 2020, showProgress = FALSE)) |>
  mutate(id_municipio = as.character(code_muni)) |>
  filter(id_municipio %in% dados_mapa2$id_municipio) |>
  left_join(dados_mapa2, by = "id_municipio")

pal_continua <- colorNumeric(
  palette = c(
    CORES_MAPA$amarelo_queimado,
    CORES_MAPA$marrom_terra,
    CORES_MAPA$vermelho_escuro),
  domain = municipios_sf2$desmatamento_ha,
  na.color = "transparent")

mapa_intensidade <- leaflet(municipios_sf2) |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(
    fillColor = ~pal_continua(desmatamento_ha),
    fillOpacity = 0.8,
    color = "#999999",
    weight = 0.8,
    opacity = 1,
    label = ~lapply(paste0(
      "<strong style='font-size:14px;'>", name_muni, " - ", sigla_uf, "</strong><br/>",
      "<span style='font-size:12px;'>",
      "Desmatamento: <strong>", format(round(desmatamento_ha, 0), big.mark = "."), " ha</strong><br/>",
      "Produção: ", format(round(producao_total_ton / 1000, 1), big.mark = "."), " mil ton<br/>",
      "Risco/ton: ", round(risco_desmat_ton, 3), " ha/ton",
      "</span>"), htmltools::HTML),
    highlightOptions = highlightOptions(
      weight = 2,
      color = "#333333",
      fillOpacity = 0.95,
      bringToFront = TRUE)) |>
  addLegend(
    position = "bottomright",
    pal = pal_continua,
    values = ~desmatamento_ha,
    title = "Desmatamento (ha)",
    opacity = 1,
    labFormat = labelFormat(big.mark = ".")) |>
  addControl(
    html = paste0(
      "<div style='background:white; padding:15px; border-radius:8px; ",
      "box-shadow: 0 2px 10px rgba(0,0,0,0.2); font-family:Arial, sans-serif;'>",
      "<h3 style='margin:0 0 10px 0; color:#333; font-size:16px;'>",
      "Intensidade de Desmatamento</h3>",
      "<p style='margin:0; color:#666; font-size:13px;'>",
      "Municípios acima do percentil ", percentil_min * 100, " (", ano_ref, ")</p>",
      "</div>"),
    position = "topleft")

htmlwidgets::saveWidget(
  mapa_intensidade,
  "outputs/maps/mapa_2_intensidade.html",
  selfcontained = TRUE)


# Mapa 2.1 estático

# Parâmetros 
ano_ref       <- 2021
percentil_min <- 0.70

# Agregação dos dados
dados_mapa2 <- df_completo |>
  filter(ano == ano_ref) |>
  group_by(id_municipio, sigla_uf) |>
  summarise(
    desmatamento_ha     = sum(desmatamento_ha, na.rm = TRUE),
    producao_total_ton  = sum(producao_total_ton, na.rm = TRUE),
    risco_desmat_ton    = mean(risco_desmat_ton, na.rm = TRUE),
    .groups = "drop")

# Percentil de corte
limiar_desmat <- quantile(dados_mapa2$desmatamento_ha, percentil_min, na.rm = TRUE)

# Base espacial – TODOS os municípios 
municipios_sf2 <- suppressMessages(
  geobr::read_municipality(year = 2020, showProgress = FALSE)) |>
  mutate(id_municipio = as.character(code_muni)) |>
  left_join(dados_mapa2, by = "id_municipio") |>
  mutate(
    destaque = if_else(desmatamento_ha >= limiar_desmat, "Acima do percentil", "Demais"))

# 5) Mapa estático -----------------------------------------------------
mapa_intensidade_static <- ggplot(municipios_sf2) +
  
  # Brasil inteiro como contexto
  geom_sf(
    aes(fill = desmatamento_ha),
    color = "white",
    linewidth = 0.08) +
  
  scale_fill_gradientn(
    colors = c(
      CORES_MAPA$amarelo_queimado,
      CORES_MAPA$marrom_terra,
      CORES_MAPA$vermelho_escuro),
    values = rescale(c(limiar_desmat, max(municipios_sf2$desmatamento_ha, na.rm = TRUE))),
    na.value = "grey90",
    name = "Desmatamento (ha)",
    labels = label_number(big.mark = ".", decimal.mark = ",")) +
  
  labs(
    title = "Intensidade de Desmatamento por Município",
    subtitle = paste0(
      "Municípios acima do percentil ",
      percentil_min * 100,
      " – ",
      ano_ref),
    caption = "Fonte: IBGE (PAM) e INPE (PRODES)") +
  
  coord_sf(expand = FALSE) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank())

# Plot
mapa_intensidade_static


# Mapa 3: Eficiência produtiva (risco/produção) ---------------------------

prod_min <- 1000
n_por_categoria <- 30

dados_eficiencia <- df_completo |>
  filter(
    ano == ano_ref,
    producao_total_ton > prod_min,
    desmatamento_ha > 0) |>
  group_by(id_municipio, sigla_uf) |>
  summarise(
    desmatamento_ha   = sum(desmatamento_ha, na.rm = TRUE),
    producao_total_ton = sum(producao_total_ton, na.rm = TRUE),
    risco_desmat_ton = mean(risco_desmat_ton, na.rm = TRUE),
    .groups = "drop") |>
  mutate(
    categoria_eficiencia = case_when(
      risco_desmat_ton < quantile(risco_desmat_ton, 0.25, na.rm = TRUE) ~ "Muito Eficiente",
      risco_desmat_ton < quantile(risco_desmat_ton, 0.50, na.rm = TRUE) ~ "Eficiente",
      risco_desmat_ton < quantile(risco_desmat_ton, 0.75, na.rm = TRUE) ~ "Ineficiente",
      TRUE ~ "Muito Ineficiente"),
    categoria_eficiencia = factor(
      categoria_eficiencia,
      levels = c("Muito Eficiente", "Eficiente", "Ineficiente", "Muito Ineficiente")))

set.seed(123)
amostra_municipios <- dados_eficiencia |>
  group_by(categoria_eficiencia) |>
  slice_max(order_by = producao_total_ton, n = n_por_categoria) |>
  ungroup()

municipios_sf3 <- suppressMessages(
  geobr::read_municipality(year = 2020, showProgress = FALSE)) |>
  mutate(id_municipio = as.character(code_muni)) |>
  filter(id_municipio %in% amostra_municipios$id_municipio) |>
  left_join(amostra_municipios, by = "id_municipio")

pal_eficiencia <- colorFactor(
  palette = c(
    CORES_MAPA$verde_floresta,
    CORES_MAPA$verde_seco,
    CORES_MAPA$amarelo_queimado,
    CORES_MAPA$vermelho_escuro),
  domain = c("Muito Eficiente", "Eficiente", "Ineficiente", "Muito Ineficiente"))

mapa_eficiencia <- leaflet(municipios_sf3) |>
  addProviderTiles(providers$Esri.WorldGrayCanvas) |>
  addPolygons(
    fillColor = ~pal_eficiencia(categoria_eficiencia),
    fillOpacity = 0.85,
    color = "#FFFFFF",
    weight = 1.5,
    opacity = 1,
    label = ~lapply(paste0(
      "<strong style='font-size:14px;'>", name_muni, " - ", sigla_uf, "</strong><br/>",
      "<span style='font-size:12px;'>",
      "Categoria: <strong>", categoria_eficiencia, "</strong><br/>",
      "Risco: ", round(risco_desmat_ton, 3), " ha/ton<br/>",
      "Desmatamento: ", format(round(desmatamento_ha, 0), big.mark = "."), " ha<br/>",
      "Produção: ", format(round(producao_total_ton / 1000, 1), big.mark = "."), " mil ton",
      "</span>"), htmltools::HTML),
    highlightOptions = highlightOptions(
      weight = 3,
      color = "#000000",
      fillOpacity = 0.95,
      bringToFront = TRUE)) |>
  addLegend(
    position = "bottomright",
    pal = pal_eficiencia,
    values = ~categoria_eficiencia,
    title = "Eficiência Produtiva",
    opacity = 1) |>
  addControl(
    html = paste0(
      "<div style='background:white; padding:15px; border-radius:8px; ",
      "box-shadow: 0 2px 10px rgba(0,0,0,0.2); font-family:Arial, sans-serif;'>",
      "<h3 style='margin:0 0 10px 0; color:#333; font-size:16px;'>",
      "Eficiência Produtiva</h3>",
      "<p style='margin:0 0 8px 0; color:#666; font-size:13px;'>",
      "Relação Desmatamento/Produção (", ano_ref, ")</p>",
      "<p style='margin:0; color:#888; font-size:11px;'>",
      "Verde = Baixo desmatamento/tonelada<br/>",
      "Vermelho = Alto desmatamento/tonelada</p>",
      "</div>"),
    position = "topleft")

htmlwidgets::saveWidget(
  mapa_eficiencia,
  "outputs/maps/mapa_3_eficiencia.html",
  selfcontained = TRUE)



# Mapa 3.1 estático

# Parâmetros 
ano_ref          <- 2021
prod_min         <- 1000
n_por_categoria  <- 30

# Cálculo da eficiência
dados_eficiencia <- df_completo |>
  filter(
    ano == ano_ref,
    producao_total_ton > prod_min,
    desmatamento_ha > 0) |>
  group_by(id_municipio, sigla_uf) |>
  summarise(
    desmatamento_ha    = sum(desmatamento_ha, na.rm = TRUE),
    producao_total_ton = sum(producao_total_ton, na.rm = TRUE),
    risco_desmat_ton   = mean(risco_desmat_ton, na.rm = TRUE),
    .groups = "drop") |>
  mutate(
    categoria_eficiencia = case_when(
      risco_desmat_ton < quantile(risco_desmat_ton, 0.25, na.rm = TRUE) ~ "Muito Eficiente",
      risco_desmat_ton < quantile(risco_desmat_ton, 0.50, na.rm = TRUE) ~ "Eficiente",
      risco_desmat_ton < quantile(risco_desmat_ton, 0.75, na.rm = TRUE) ~ "Ineficiente",
      TRUE ~ "Muito Ineficiente"),
    categoria_eficiencia = factor(
      categoria_eficiencia,
      levels = c(
        "Muito Eficiente",
        "Eficiente",
        "Ineficiente",
        "Muito Ineficiente")))

# Amostra por categoria 
set.seed(123)

amostra_municipios <- dados_eficiencia |>
  group_by(categoria_eficiencia) |>
  slice_max(order_by = producao_total_ton, n = n_por_categoria) |>
  ungroup() |>
  mutate(id_municipio = as.character(id_municipio))

# Base espacial – TODOS os municípios
municipios_sf3 <- suppressMessages(
  geobr::read_municipality(year = 2020, showProgress = FALSE)) |>
  mutate(id_municipio = as.character(code_muni)) |>
  left_join(amostra_municipios, by = "id_municipio")

# Mapa estático
mapa_eficiencia_static <- ggplot(municipios_sf3) +
  
  geom_sf(
    aes(fill = categoria_eficiencia),
    color = "white",
    linewidth = 0.1) +
  
  scale_fill_manual(
    values = c(
      "Muito Eficiente"   = CORES_MAPA$verde_floresta,
      "Eficiente"         = CORES_MAPA$verde_seco,
      "Ineficiente"       = CORES_MAPA$amarelo_queimado,
      "Muito Ineficiente" = CORES_MAPA$vermelho_escuro),
    na.value = "grey90",
    name = "Eficiência produtiva") +
  
  labs(
    title = "Eficiência Produtiva dos Municípios",
    subtitle = paste0(
      "Relação desmatamento / produção – ",
      ano_ref,
      " | Amostra: ",
      n_por_categoria,
      " municípios por categoria"),
    caption = "Fonte: IBGE (PAM) e INPE (PRODES)") +
  
  coord_sf(expand = FALSE) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank())

# Plot 
mapa_eficiencia_static


# Mapa 4: Análise por Estado ----------------------------------------------

dados_uf <- df_completo |>
  filter(ano == ano_ref) |>
  group_by(sigla_uf) |>
  summarise(
    producao_total     = sum(producao_total_ton, na.rm = TRUE),
    desmatamento_total = sum(desmatamento_ha, na.rm = TRUE),
    n_municipios       = n_distinct(id_municipio),
    risco_medio        = mean(risco_desmat_ton, na.rm = TRUE),
    .groups = "drop") |>
  mutate(
    producao_milhoes_ton = producao_total / 1000000,
    intensidade_desmat = desmatamento_total / (producao_total / 1000))

estados_sf <- suppressMessages(
  geobr::read_state(year = 2020, showProgress = FALSE)) |>
  left_join(dados_uf, by = c("abbrev_state" = "sigla_uf"))

pal_intensidade <- colorNumeric(
  palette = c(
    CORES_MAPA$verde_floresta,
    CORES_MAPA$amarelo_queimado,
    CORES_MAPA$marrom_terra,
    CORES_MAPA$vermelho_escuro),
  domain = estados_sf$intensidade_desmat,
  na.color = "#EEEEEE")

centroides <- estados_sf |>
  st_centroid() |>
  filter(!is.na(producao_total))

mapa_estados <- leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(
    data = estados_sf,
    fillColor = ~pal_intensidade(intensidade_desmat),
    fillOpacity = 0.75,
    color = "#666666",
    weight = 1.5,
    label = ~lapply(paste0(
      "<strong style='font-size:14px;'>", name_state, "</strong><br/>",
      "<span style='font-size:12px;'>",
      "Produção: ", round(producao_milhoes_ton, 2), " milhões ton<br/>",
      "Desmatamento: ", format(round(desmatamento_total, 0), big.mark = "."), " ha<br/>",
      "Intensidade: ", round(intensidade_desmat, 2), " ha/mil ton<br/>",
      "Municípios: ", n_municipios,
      "</span>"), htmltools::HTML)) |>
  addCircleMarkers(
    data = centroides,
    radius = ~sqrt(producao_total / 100000) * 5,
    fillColor = CORES_GRAFICO$azul_medio,
    fillOpacity = 0.7,
    color = CORES_GRAFICO$azul_claro,
    weight = 2,
    label = ~lapply(paste0(
      "<strong style='font-size:14px;'>", name_state, "</strong><br/>",
      "<span style='font-size:12px;'>",
      "Produção: ", round(producao_milhoes_ton, 2), " milhões ton",
      "</span>"), htmltools::HTML)) |>
  addLegend(
    position = "bottomright",
    pal = pal_intensidade,
    values = estados_sf$intensidade_desmat,
    title = "Intensidade<br/>(ha/mil ton)",
    labFormat = labelFormat(digits = 1)) |>
  addControl(
    html = paste0(
      "<div style='background:white; padding:15px; border-radius:8px; ",
      "box-shadow: 0 2px 10px rgba(0,0,0,0.2); font-family:Arial, sans-serif;'>",
      "<h3 style='margin:0 0 10px 0; color:#333; font-size:16px;'>",
      "Produção e Desmatamento por Estado (", ano_ref, ")</h3>",
      "<p style='margin:0; color:#666; font-size:12px;'>",
      "Cor: Intensidade de desmatamento<br/>",
      "Círculos: Volume de produção</p>",
      "</div>"),
    position = "topleft")

htmlwidgets::saveWidget(
  mapa_estados,
  "outputs/maps/mapa_4_estados.html",
  selfcontained = TRUE)


# Mapa 4.1 estático

# Parâmetros
ano_ref <- 2021

# Agregação por UF 
dados_uf <- df_completo |>
  filter(ano == ano_ref) |>
  group_by(sigla_uf) |>
  summarise(
    producao_total     = sum(producao_total_ton, na.rm = TRUE),
    desmatamento_total = sum(desmatamento_ha, na.rm = TRUE),
    n_municipios       = n_distinct(id_municipio),
    risco_medio        = mean(risco_desmat_ton, na.rm = TRUE),
    .groups = "drop") |>
  mutate(
    producao_milhoes_ton = producao_total / 1e6,
    intensidade_desmat  = desmatamento_total / (producao_total / 1000))

# Base espacial – Estados
estados_sf <- suppressMessages(
  geobr::read_state(year = 2020, showProgress = FALSE)) |>
  left_join(dados_uf, by = c("abbrev_state" = "sigla_uf"))

# Centroides para pontos de produção
centroides <- estados_sf |>
  st_centroid() |>
  filter(!is.na(producao_total))

# Mapa estático
mapa_estados_static <- ggplot(estados_sf) +
  
  # Estados coloridos pela intensidade
  geom_sf(
    aes(fill = intensidade_desmat),
    color = "white",
    linewidth = 0.2) +
  
  # Pontos proporcionais à produção
  geom_sf(
    data = centroides,
    aes(size = producao_milhoes_ton),
    color = CORES_GRAFICO$azul_claro,
    fill  = CORES_GRAFICO$azul_medio,
    shape = 21,
    alpha = 0.7,
    stroke = 0.8) +
  
  scale_fill_gradientn(
    colors = c(
      CORES_MAPA$verde_floresta,
      CORES_MAPA$amarelo_queimado,
      CORES_MAPA$marrom_terra,
      CORES_MAPA$vermelho_escuro),
    na.value = "grey90",
    name = "Intensidade de desmatamento\n(ha / mil ton)",
    labels = label_number(decimal.mark = ",", big.mark = ".")) +
  
  scale_size_continuous(
    range = c(3, 12),
    name = "Produção\n(milhões ton)") +
  
  labs(
    title = "Produção e Desmatamento por Estado",
    subtitle = paste0(
      "Cor: intensidade de desmatamento | Tamanho: volume de produção – ",
      ano_ref
    ),
    caption = "Fonte: IBGE (PAM) e INPE (PRODES)") +
  
  coord_sf(expand = FALSE) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank())

# Plot
mapa_estados_static

# 8. Exportações dos resultados -------------------------------------------


# Criar diretórios de saída para tabelas
if (!dir.exists("outputs/tables")) dir.create("outputs/tables")

# Dataset completo processado
readr::write_csv(
  df_completo,
  "outputs/tables/dataset_integrado.csv")

# Predições do modelo no conjunto de teste
readr::write_csv(
  avaliacao_modelo$resultados,
  "outputs/tables/predicoes_teste.csv")

# 50 municípios com maior desmatamento em 2021
top50 <- df_completo |>
  filter(ano == 2021) |>
  arrange(desc(desmatamento_ha)) |>
  slice_head(n = 50) |>
  select(
    id_municipio, sigla_uf, ano,
    desmatamento_ha, producao_total_ton, area_colhida_ha,
    produtividade, classe_risco, risco_desmat_ton)

readr::write_csv(
  top50,
  "outputs/tables/top50_municipios_2021.csv")

# Análises descritivas
readr::write_csv(
  analise_descritiva$por_classe_risco,
  "outputs/tables/analise_descritiva_por_risco.csv")

readr::write_csv(
  analise_descritiva$temporal,
  "outputs/tables/analise_descritiva_temporal.csv")

# Métricas do modelo
readr::write_csv(
  metricas_complementares,
  "outputs/tables/metricas_modelo.csv")

# 9. Reporte final --------------------------------------------------------

# Criar diretórios de saída para tabelas
if (!dir.exists("outputs/tables")) dir.create("outputs/tables")

relatorio_final <- tibble(
  Item = c(
    "Observações no dataset",
    "Municípios únicos",
    "Período analisado",
    "Features utilizadas",
    "MAE (hectares)",
    "RMSE (hectares)",
    "R² do modelo",
    "MAPE (%)"),
  Valor = c(
    format(nrow(df_completo), big.mark = "."),
    format(n_distinct(df_completo$id_municipio), big.mark = "."),
    sprintf("%d-%d", min(df_completo$ano), max(df_completo$ano)),
    length(FEATURES),
    format(metricas_complementares$MAE, big.mark = ".", decimal.mark = ","),
    format(metricas_complementares$RMSE, big.mark = ".", decimal.mark = ","),
    metricas_complementares$R2,
    format(metricas_complementares$MAPE, big.mark = ".", decimal.mark = ",")))

readr::write_csv(
  relatorio_final,
  "outputs/tables/relatorio_execucao.csv")
view(relatorio_final)


# Mapas estáticos ---------------------------------------------------------

# Diretório
dir.create("outputs/maps_estaticos", recursive = TRUE, showWarnings = FALSE)

# Parâmetros
largura  <- 12   # polegadas
altura   <- 8
dpi_plot <- 300  # padrão artigo / relatório

# 1) Mapa 1.1 – Top municípios
ggsave(
  filename = "outputs/maps_estaticos/mapa_1_top_municipios_desmatamento.png",
  plot     = mapa_top_risco_static,
  width    = largura,
  height   = altura,
  dpi      = dpi_plot,
  bg       = "white")

# 2) Mapa 2.1 – Intensidade de desmatamento
ggsave(
  filename = "outputs/maps_estaticos/mapa_2_intensidade_desmatamento.png",
  plot     = mapa_intensidade_static,
  width    = largura,
  height   = altura,
  dpi      = dpi_plot,
  bg       = "white")

# 3) Mapa 3.1 – Eficiência produtiva
ggsave(
  filename = "outputs/maps_estaticos/mapa_3_eficiencia_produtiva.png",
  plot     = mapa_eficiencia_static,
  width    = largura,
  height   = altura,
  dpi      = dpi_plot,
  bg       = "white")

# 4) Mapa 4.1 – Produção e desmatamento por estado -----------------------
ggsave(
  filename = "outputs/maps_estaticos/mapa_4_producao_desmatamento_estados.png",
  plot     = mapa_estados_static,
  width    = largura,
  height   = altura,
  dpi      = dpi_plot,
  bg       = "white")