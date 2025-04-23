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

# Malha de regiões de Campinas:
malha_campinas <- read_sf("Malha/Campinas2.shp") %>%
    select(APG, POP_2022)

# Função para carregamento e tratamento das bases de dados:
func_carregamento <- function(string) {
    df <- read_xlsx(string) %>%
        select(NOME_DELEGACIA, NUM_BO, NOME_MUNICIPIO,
               LATITUDE, LONGITUDE, ANO) %>%
        filter(NOME_MUNICIPIO == "CAMPINAS") %>%
        mutate(ID = str_c(NOME_DELEGACIA, NUM_BO)) %>%
        distinct(ID, .keep_all = T) %>%
        mutate(LATITUDE = as.numeric(LATITUDE),
               LONGITUDE = as.numeric(LONGITUDE)) %>%
        filter(LATITUDE != 0)
}

# Nomes das bases para função:
nomes_bases <- paste0("~Downloads/Bases/Bayesiana/", c("CelularesSubtraidos_2017.xlsx", "CelularesSubtraidos_2018.xlsx",
                                                       "CelularesSubtraidos_2019.xlsx", "CelularesSubtraidos_2020.xlsx",
                                                       "CelularesSubtraidos_2021.xlsx", "CelularesSubtraidos_2022.xlsx",
                                                       "CelularesSubtraidos_2023.xlsx", "CelularesSubtraidos_2024.xlsx"))

# Lista com roubos em Campinas, por ano:
lista_roubos <- map(nomes_bases, func_carregamento)

# Função para extração e transformação das coordenadas:
func_lat_long <- function(df) {
    out <- df %>%
        select(ANO, LATITUDE, LONGITUDE) %>%
        st_as_sf(coords = c("LONGITUDE", "LATITUDE"),
                 crs = 4326)
    out <- st_transform(out, crs = 31983)
}

# Posições pontuais dos roubos:
locais_roubos <- map(lista_roubos, func_lat_long)

# Função para contagem da interseção espacial entre a malha e os pontos:
func_contagem <- function(df, malha) {
    ano <- df$ANO[1]
    out <- malha %>%
        mutate(Roubos = lengths(st_intersects(malha, df)),
               Ano = ano)
    return(out)
}

# Aplicação da função para todos os anos:
contagem_campinas <- map(locais_roubos, func_contagem, malha = malha_campinas)

# Transformação em um data.frame unificado:
contagem_campinas <- do.call(rbind, contagem_campinas)

# Gráfico com a evolução de roubos ao longo dos anos:
contagem_campinas %>%
    ggplot(aes(fill = Roubos)) +
    geom_sf() +
    geom_sf_label(aes(label = Roubos), size = 1) +
    scale_fill_distiller(palette = "Spectral") +
    facet_wrap(~Ano)

#----Modelagem----

# Função para calcular a expectativa de roubos:
func_expect <- function(df){
    expected(df$POP_2022, df$Roubos, 1)
}

E_campinas <- numeric()

for (i in 2017:2024) {
    E_campinas <- contagem_campinas %>%
        filter(Ano == i) %>%
        func_expect() %>%
        c(E_campinas, .)
}

# Risco Relativo a priori:
contagem_campinas$RRAPriori <- contagem_campinas$Roubos/E_campinas

contagem_campinas %>%
    ggplot(aes(fill = RRAPriori)) +
    geom_sf() +
    scale_fill_distiller(palette = "Spectral") +
    facet_wrap(~Ano)

# Matriz de adjacência do modelo:
adj_campinas <- poly2nb(malha_campinas)
mat_adj <- as(nb2mat(adj_campinas, style = "B"), "Matrix")

# Mudança dos anos:
contagem_campinas$Ano2 <- contagem_campinas$Ano - 2017 + 1

# Efeitos aleatórios espaço-temporais:
contagem_campinas$ea_u <- rep(1:17, 8)

# Modelo:
modelo_completo <- inla(Roubos ~ 1 +
                            f(ea_u, model = "bym", graph = mat_adj,
                              group = Ano2, control.group = list(model = "ar1")),
                        data = contagem_campinas, E = E_campinas, family = "nbinomial",
                        control.predictor = list(compute = TRUE, link = 1))

# Valores ajustados para cada bairro e ano:
modelo_completo$summary.fitted.values

# Média do Risco Relativo a Posteriori:
contagem_campinas$RRAP <- modelo_completo$summary.fitted.values[, "mean"]

# Distribuição espaço-temporal do risco relativo a posteriori:
contagem_campinas %>%
    ggplot(aes(fill = RRAP)) +
    geom_sf() +
    scale_fill_distiller(palette = "Spectral") +
    facet_wrap(~Ano)

# Diferença no risco relativo a priori e posteriori para Barão Geraldo:
contagem_campinas %>%
    filter(APG == "Barão Geraldo") %>%
    select(Ano, RRAP, RRAPriori) %>%
    ggplot(aes(x = Ano, y = RRAP, colour = "Priori")) +
    geom_line() +
    geom_point() +
    geom_point(aes(y = RRAPriori, colour = "Posteriori")) +
    geom_line(aes(y = RRAPriori, colour = "Posteriori")) +
    scale_colour_manual(values = c("darkgreen", "tomato")) +
    labs(x = "Anos", y = "Risco Relativo")

#----Mapa Animado----

# Pacote para criação de animações com ggplot:
library(gganimate)

# O seguinte código cria 100 imagens no diretório:
ggplot(contagem_campinas) +
    geom_sf(aes(fill = RAP)) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()) +
    scale_fill_distiller(palette = "Spectral") +
    transition_time(Ano) +
    labs(title = "Ano: {round(frame_time, 0)}")
