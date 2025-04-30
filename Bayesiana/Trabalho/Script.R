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

# Função para carregamento e tratamento das bases de dados:
func_carregamento_sp <- function(string) {
    df <- read_xlsx(string) %>%
        select(NOME_DELEGACIA, NUM_BO, NOME_MUNICIPIO,
               LATITUDE, LONGITUDE, ANO) %>%
        mutate(ID = str_c(NOME_DELEGACIA, NUM_BO)) %>%
        distinct(ID, .keep_all = T) %>%
        group_by(NOME_MUNICIPIO) %>%
        summarise(Total = n(),
                  Ano = mean(ANO))
}

# Nomes das bases para função:
nomes_bases <- paste0("~/Downloads/Bases/Bayesiana/", c("CelularesSubtraidos_2017.xlsx", "CelularesSubtraidos_2018.xlsx",
                                                        "CelularesSubtraidos_2019.xlsx", "CelularesSubtraidos_2020.xlsx",
                                                        "CelularesSubtraidos_2021.xlsx", "CelularesSubtraidos_2022.xlsx",
                                                        "CelularesSubtraidos_2023.xlsx", "CelularesSubtraidos_2024.xlsx"))

# Lista com 8 bases de dados:
bases_sp <- map(nomes_bases, func_carregamento_sp)

# Pivoteamento para ano:
roubos_sp <- do.call(rbind, bases_sp) %>%
    pivot_wider(names_from = Ano, values_from = Total, values_fill = 0)

# Criar arquivo para manipulação manual dos nomes dos municípios:
roubos_sp %>%
    arrange(NOME_MUNICIPIO) %>%
    write_csv(file = "roubos_sp.csv")

# Arquivo ajustado:
roubos_sp_total <- read.csv("roubos_sp_total.csv") %>%
    select(NM_MUN = NOME_MUNICIPIO, everything())

# Malha municipal:
mapa_sp_total <- read_sf("Malha/SP_Municipios_2023.shp")

#----Modelo espaço-temporal----

# Matriz de adjacência entre os municípios:
nb_sp <- poly2nb(mapa_sp_total)

nb2INLA("mapa_sp.adj", nb_sp)

# Grafo para INLA:
g <- inla.read.graph(filename = "mapa_sp.adj")

# Conexão das bases com o mapa:
mapa_sp_total <- mapa_sp_total %>%
    left_join(municipios_sp %>% select(-NM_MUN), by = "CD_MUN")

mapa_sp_total <- mapa_sp_total %>%
    select(NM_MUN, Populacao, Escolarizacao, IDHM) %>%
    left_join(roubos_sp_total, by = "NM_MUN") %>%
    pivot_longer(cols = !c(NM_MUN, geometry, Populacao, Escolarizacao, IDHM), names_to = "Ano", values_to = "Roubos") %>%
    mutate(Ano = as.numeric(str_sub(Ano, start = 2))) %>%
    arrange(Ano, NM_MUN)

# Valores esperados de assaltos por ano:
E_sp_total <- numeric()

for (i in 2017:2024) {
    E_sp_total <- expected(
    (filter(mapa_sp_total, Ano == i))$Populacao,
    (filter(mapa_sp_total, Ano == i))$Roubos,
    1
    ) %>%
        c(E_sp_total, .)
}

mapa_sp_total$E <- E_sp_total

# Efeitos aleatórios:

## Espacialmente correlacionados:
mapa_sp_total$ea_u <- rep(1:nrow(roubos_sp_total), 8)

## Temporais
mapa_sp_total$Ano2 <- mapa_sp_total$Ano - 2017 + 1

# Adição de 2025 para predição:
sp_2025 <- data.frame(
    NM_MUN = mapa_sp_total$NM_MUN[1:nrow(roubos_sp_total)],
    Populacao = mapa_sp_total$Populacao[1:nrow(roubos_sp_total)],
    Escolarizacao = mapa_sp_total$Escolarizacao[1:nrow(roubos_sp_total)],
    IDHM = mapa_sp_total$IDHM[1:nrow(roubos_sp_total)],
    geometry = mapa_sp_total$geometry[1:nrow(roubos_sp_total)],
    Ano = 2025,
    Roubos = NA,
    E = NA,
    ea_u = 1:nrow(roubos_sp_total),
    Ano2 = 9
) %>%
    st_as_sf() %>%
    st_transform(crs = 4674)

mapa_sp_total <- rbind(mapa_sp_total, sp_2025)

#----Modelagem----

# Fórmula usando modelo AR1 para temporal e Besag-York-Mollier combinado
# Com ano AR1 para espaço-temporal:
formula <- Roubos ~ IDHM + Escolarizacao +
    f(Ano, model = "ar1") +
    f(ea_u, model = "bym", graph = g, scale.model = T,
      group = Ano2, control.group = list(model = "ar1"))

# Modelo final:
modelo_sp_total <- inla(formula, family = "poisson", data = mapa_sp_total, E = E,
                        control.predictor = list(compute = T, link = 1),
                        control.compute = list(return.marginals.predictor = T,
                                               dic = T))

# Estimativas pontuais dos betas:
modelo_sp_total$summary.fixed

# Distribuições a posteriori dos betas:
lista_betas_sp <- modelo_sp_total$marginals.fixed

posterioris_betas_sp <- vector(mode = "list", length = length(lista_betas_sp))

for (i in 1:length(lista_betas_sp)) {
    posterioris_betas_sp[[i]] <- lista_betas_sp[[i]] %>%
        ggplot(aes(x = x, y = y)) +
        geom_line() +
        labs(title = names(lista_betas_sp)[i])
}

ggpubr::ggarrange(plotlist = posterioris_betas_sp)

# Distribuições a posteriori para os valores dos hiper-parâmetros:
lista_hiperparametros_sp <- modelo_sp_total$marginals.hyperpar

posterioris_hiper_sp <- vector(mode = "list", length = length(lista_hiperparametros_sp))

for (i in 1:length(lista_hiperparametros_sp)) {
    posterioris_hiper_sp[[i]] <- lista_hiperparametros_sp[[i]] %>%
        ggplot(aes(x = x, y = y)) +
        geom_line() +
        labs(title = names(lista_hiperparametros_sp)[i])
}

ggpubr::ggarrange(plotlist = posterioris_hiper_sp)

# Valores ajustados para cada cidade e ano:
modelo_sp_total$summary.fitted.values

# Média do Risco Relativo a Posteriori:
mapa_sp_total$RRAP <- modelo_sp_total$summary.fitted.values[, "mean"]

mapa_sp_total$RRAP_L <- modelo_sp_total$summary.fitted.values[, "0.025quant"]
mapa_sp_total$RRAP_U <- modelo_sp_total$summary.fitted.values[, "0.975quant"]

# Distribuição espaço-temporal do risco relativo a posteriori:
mapa_sp_total %>%
    ggplot(aes(fill = RRAP)) +
    geom_sf() +
    scale_fill_distiller(palette = "Spectral",
                         limits = c(min(mapa_sp_total$RRAP_L),
                                    max(mapa_sp_total$RRAP_U))) +
    facet_wrap(~Ano) +
    ggpubr::theme_classic2()

mapa_sp_total %>%
    ggplot(aes(fill = RRAP_L)) +
    geom_sf() +
    scale_fill_distiller(palette = "Spectral",
                         limits = c(min(mapa_sp_total$RRAP_L),
                                    max(mapa_sp_total$RRAP_U))) +
    facet_wrap(~Ano) +
    ggpubr::theme_classic2()

mapa_sp_total %>%
    ggplot(aes(fill = RRAP_U)) +
    geom_sf() +
    scale_fill_distiller(palette = "Spectral",
                         limits = c(min(mapa_sp_total$RRAP_L),
                                    max(mapa_sp_total$RRAP_U))) +
    facet_wrap(~Ano) +
    ggpubr::theme_classic2()

# Predição a Posteriori:

# Densidades marginais do risco relativo a posteriori para a cidade de São Paulo:
marginais_sp <- list()

for (i in 1:8) {
    marginais_sp[[i]] <- inla.smarginal(modelo_sp_total$marginals.fitted.values[[573 + nrow(roubos_sp_total) * (i - 1)]]) %>%
        data.frame() %>%
        mutate(Ano = 2017 + (i - 1))
}

marginais_sp <- do.call(rbind, marginais_sp)

# Mudança no risco relativo a posteriori para São Paulo:
marginais_sp %>%
    ggplot(aes(x = x, y = y)) +
    geom_line() +
    facet_wrap(~Ano) +
    labs(x = "Risco a posteriori", y = "Densidade") +
    theme_bw()
