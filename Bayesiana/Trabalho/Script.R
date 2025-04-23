#----Pacotes exigidos----

# Manipulação dos dados:
library(tidyverse)

# Manipulação espacial:
library(sf)

# Criação de matriz de adjacência de polígonos:
library(spdep)

# Cálculo de quantidade esperada empírica de roubos por município:
library(SpatialEpi)

# Carregamento dos dados:
library(readxl)

# Modelagem:
library(INLA)

#----Manipulação dos dados----

# Carregamento dos dados:
celulares_2023 <- read_xlsx("CelularesSubtraidos_2023.xlsx") %>%
    select(NOME_DELEGACIA, NUM_BO, NOME_MUNICIPIO,
           LATITUDE, LONGITUDE, MES)

# Criar coluna para identificação de observações duplicadas,
# conforme a metodologia no dicionário dos dados:
celulares_2023 <- celulares_2023 %>%
    mutate(ID = str_c(NOME_DELEGACIA, NUM_BO))

# Manter as observações únicas:
celulares_2023 <- celulares_2023 %>%
    distinct(ID, .keep_all = T)

# Informações municipais:
municipios_sp <- read_xlsx("municipios_sp.xlsx") %>%
    mutate(CD_MUN = as.character(CD_MUN))

# Cria arquivo com os nomes dos municípios e o total de roubos:
celulares_2023 %>%
    group_by(NOME_MUNICIPIO) %>%
    summarise(Total = n()) %>%
    arrange(NOME_MUNICIPIO) %>%
    write_csv(file = "roubos_por_municipio.csv")

# Carrega arquivo com os nomes ajustados manualmente:
roubos_2023 <- read.csv("roubos_total.csv")

# Malha municipal:
mapa_sp <- read_sf("Malha/SP_Municipios_2023.shp")

mapa_sp <- mapa_sp %>%
    left_join(municipios_sp %>% select(-NM_MUN), by = "CD_MUN")

# Junção dos dados no mapa:
mapa_sp <- mapa_sp %>%
    left_join(roubos_2023, by = "NM_MUN") %>%
    mutate(ROUBOS = case_when(is.na(ROUBOS) ~ 0,
                              .default = ROUBOS))

E <- expected(mapa_sp$Populacao, mapa_sp$ROUBOS, 1)

mapa_sp$E <- E

nb <- poly2nb(mapa_sp)

nb2INLA("mapa_sp.adj", nb)

g <- inla.read.graph(filename = "mapa_sp.adj")

# Efeitos aleatórios:
## Espacialmente correlacionados:
mapa_sp$ea_u <- 1:nrow(mapa_sp)

## Não correlacionados:
mapa_sp$ea_v <- 1:nrow(mapa_sp)

#----Modelagem----

formula <- ROUBOS ~ IDHM + Escolarizacao +
    f(ea_u, model = "besag", graph = g, scale.model = T) +
    f(ea_v, model = "iid")

modelo <- inla(formula, family = "poisson", data = mapa_sp, E = E,
               control.predictor = list(compute = T),
               control.compute = list(return.marginals.predictor = T))

modelo$summary.fixed

# Estimação do risco relativo a posteriori:
mapa_sp$RAP <- modelo$summary.fitted.values[, "mean"]

modelo$summary.fitted.values

mapa_sp %>%
    select(RAP) %>%
    plot()

#----Testes----

# Malha de bairros em Campinas:
malhas <- read_sf("Malha/Campinas2.shp") %>%
    select(APG, POP_2022)

plot(malhas)

# Alterando para valores numericos as coordenadas:
celulares_2023$LATITUDE <- as.numeric(celulares_2023$LATITUDE, digits = 10)
celulares_2023$LONGITUDE <- as.numeric(celulares_2023$LONGITUDE, digits = 10)

campinas_2023 <- subset(celulares_2023, NOME_MUNICIPIO == "CAMPINAS" & LATITUDE != 0)

roubos_campinas_2023 <- campinas_2023 %>%
    dplyr::select(LATITUDE, LONGITUDE) %>%
    st_as_sf(coords = c("LONGITUDE", "LATITUDE"),
             crs = 4326)

roubos_campinas_2023 <- st_transform(roubos_campinas_2023, crs = 31983)

contagem_2023 <- malhas %>%
    mutate(roubos = lengths(st_intersects(malhas, roubos_campinas_2023)),
           Ano = 2023)

plot(contagem_2023["roubos"])

nb_campinas <- poly2nb(contagem_2023)

nb2INLA("mapa_campinas.adj", nb_campinas)

g_campinas <- inla.read.graph(filename = "mapa_campinas.adj")

# Efeitos aleatórios:
## Espacialmente correlacionados:
contagem_2023$ea_u <- 1:nrow(contagem)

## Não correlacionados:
contagem_2023$ea_v <- 1:nrow(contagem)

#----Modelagem----

E_campinas_2023 <- expected(contagem_2023$POP_2022, contagem_2023$roubos, 1)

formula_campinas <- roubos ~ 1 +
    f(ea_u, model = "besag", graph = g_campinas, scale.model = T) +
    f(ea_v, model = "iid")

modelo_campinas_2023 <- inla(formula_campinas, family = "poisson", data = contagem_2023, E = E_campinas_2023,
                             control.predictor = list(compute = T),
                             control.compute = list(return.marginals.predictor = T))

# Sumário dos efeitos fixos:
modelo_campinas_2023$summary.fitted.values

# Sumário dos efeitos aleatórios:
modelo_campinas_2023$summary.random

# Estimação do risco relativo a posteriori:
contagem_2023$RAP <- modelo_campinas_2023$summary.fitted.values[, "mean"]

# Mapa do Risco a Posteriori:
contagem_2023 %>%
    dplyr::select(RAP) %>%
    plot()

plot(contagem_2023[c("roubos", "RAP")])

contagem_2023 %>%
    ggplot(aes(fill = RAP)) +
    geom_sf() +
    geom_sf_label(aes(label = round(RAP, 1))) +
    scale_fill_continuous(type = "viridis")

# Predição a Posteriori:

# Densidade marginal do risco a posteriori para Barão Geraldo:
marginal_bg <- inla.smarginal(modelo_campinas$marginals.fitted.values[[3]]) %>%
    data.frame()

marginal_bg %>%
    ggplot(aes(x = x, y = y)) +
    geom_line()

# Probabilidade de Barão Geraldo ter um risco de assalto maior que 1.2:
1 - inla.pmarginal(1.2, modelo_campinas$marginals.fitted.values[[3]])

#----Espaço-temporal----

library("DClusterm")
data(brainNM)

brainst@data

celulares_2024 <- read_xlsx("CelularesSubtraidos_2024.xlsx") %>%
    dplyr::select(NOME_DELEGACIA, NUM_BO, NOME_MUNICIPIO,
                  LATITUDE, LONGITUDE, MES)

celulares_2024$LATITUDE <- as.numeric(celulares_2024$LATITUDE, digits = 10)
celulares_2024$LONGITUDE <- as.numeric(celulares_2024$LONGITUDE, digits = 10)

campinas_2024 <- subset(celulares_2024, NOME_MUNICIPIO == "CAMPINAS" & LATITUDE != 0)

roubos_campinas_2024 <- campinas_2024 %>%
    dplyr::select(LATITUDE, LONGITUDE) %>%
    st_as_sf(coords = c("LONGITUDE", "LATITUDE"),
             crs = 4326)

roubos_campinas_2024 <- st_transform(roubos_campinas_2024, crs = 31983)

contagem_2024 <- malhas %>%
    mutate(roubos = lengths(st_intersects(malhas, roubos_campinas_2024)),
           Ano = 2024)

plot(contagem_2024["roubos"])

contagem_2024$ea_u <- 1:nrow(contagem)
contagem_2024$ea_v <- 1:nrow(contagem)

E_campinas_2024 <- expected(contagem_2024$POP_2022, contagem_2024$roubos, 1)

E_completo <- c(E_campinas_2023, E_campinas_2024)

matriz_nb_campinas <- as(nb2mat(nb_campinas, style = "B"), "Matrix")

contagem_completa <- rbind(contagem_2023, contagem_2024)

formula_completa <- roubos ~ 1 + f(Ano, model = "ar1", hyper = list(prec = list(param = c(0.001, 0.001)))) +
    f(ea_u, model = "besag", graph = matriz_nb_campinas, hyper = list(prec = list(param = c(0.001, 0.001))))

modelo_completo <- inla(formula_completa, family = "poisson", data = contagem_completa, E = E_completo,
                        control.predictor = list(compute = T),
                        control.compute = list(return.marginals.predictor = T))

summary(modelo_completo)

contagem_completa$RAP <- modelo_completo$summary.fitted.values[, "mean"]

contagem_completa %>%
    ggplot(aes(fill = RAP)) +
    geom_sf() +
    geom_sf_label(aes(label = round(RAP, 1))) +
    facet_wrap(~Ano)
