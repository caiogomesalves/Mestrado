# Pacotes
library(dplyr)
library(parallel)

# Funções

dados_bras <- read.csv("~/brasileirao_2026.csv")

# Criando as taxas theta e phi:
taxas_estimadas <- left_join(
(dados_bras %>%
 mutate(time = time_mandante) %>%
 group_by(time) %>%
 summarise(theta_m = mean(gols_mandante, na.rm = T),
           phi_v = mean(gols_visitante, na.rm = T))),
(dados_bras %>%
 mutate(time = time_visitante) %>%
 group_by(time) %>%
 summarise(theta_v = mean(gols_visitante, na.rm = T),
           phi_m = mean(gols_mandante, na.rm = T))), by = "time")

cria_df_taxas <- function(df, taxas) {
    times <- unique(taxas$time)
    df_mandante <- left_join(df, (taxas %>%
                                  select(time_mandante = time,
                                         theta_m, phi_m)), by = "time_mandante") %>%
        select(theta_m, phi_m)
    df_visitante <- left_join(df, (taxas %>%
                                  select(time_visitante = time,
                                         theta_v, phi_v)), by = "time_visitante") %>%
        select(theta_v, phi_v)
    df_completo <- cbind(df, df_mandante, df_visitante)
    return(df_completo)
}

df_completo <- cria_df_taxas(dados_bras, taxas_estimadas)

# 1. Separar jogos realizados e futuros
jogos_realizados <- df_completo %>% filter(!is.na(gols_mandante))
jogos_futuros <- df_completo %>% filter(is.na(gols_mandante))

# 2. Calcular a classificação atual UMA ÚNICA VEZ considerando desempates
calcula_tabela <- function(df) {
    df %>%
        mutate(
            pts_m = case_when(gols_mandante > gols_visitante ~ 3,
                              gols_mandante == gols_visitante ~ 1,
                              TRUE ~ 0),
            pts_v = case_when(gols_visitante > gols_mandante ~ 3,
                              gols_mandante == gols_visitante ~ 1,
                              TRUE ~ 0),
            vit_m = ifelse(gols_mandante > gols_visitante, 1, 0),
            vit_v = ifelse(gols_visitante > gols_mandante, 1, 0)
        ) %>%
        pivot_longer(cols = c(time_mandante, time_visitante),
                     names_to = "mando", values_to = "time") %>%
        mutate(
            pontos = ifelse(mando == "time_mandante", pts_m, pts_v),
            vitorias = ifelse(mando == "time_mandante", vit_m, vit_v),
            gols_pro = ifelse(mando == "time_mandante", gols_mandante,
                              gols_visitante),
            gols_contra = ifelse(mando == "time_mandante", gols_visitante,
                                 gols_mandante),
            saldo = gols_pro - gols_contra
        ) %>%
        group_by(time) %>%
        summarise(across(c(pontos, vitorias, saldo, gols_pro),
                         sum, na.rm = T), .groups = "drop")
}

tabela_base <- calcula_tabela(jogos_realizados)

# 3. Função otimizada para simular apenas o restante
simula_rodadas_restantes <- function(...) {
    sim_futuros <- jogos_futuros %>%
        mutate(
            gols_mandante = rpois(n(), (theta_m + phi_v)/2),
            gols_visitante = rpois(n(), (theta_v + phi_m)/2)
        )
    tabela_simulada <- calcula_tabela(sim_futuros)
    tabela_final <- bind_rows(tabela_base, tabela_simulada) %>%
        group_by(time) %>%
        summarise(across(everything(), sum, na.rm = T), .groups = "drop") %>%
        arrange(desc(pontos), desc(vitorias), desc(saldo), desc(gols_pro)) %>%
        mutate(posicao = row_number())
    # Retornar uma lista com informações para responder a todas as perguntas
    list(
        campeao = tabela_final$time[1],
        pontos_campeao = tabela_final$pontos[1],
        rebaixados = tabela_final$time[17:20],
        pontos_todos = tabela_final %>% select(time, pontos, posicao)
    )
}

# 4. Execução paralela (tira proveito do ambiente Linux)
set.seed(5050)
resultados <- mclapply(1:100000, simula_rodadas_restantes, mc.cores = detectCores() - 1)

save(resultados, "resultados_bras.RData")
