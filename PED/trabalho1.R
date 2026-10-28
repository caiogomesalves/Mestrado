#====Pacotes necessários====

library(tidyverse)

#====Funções Auxiliares====

campeonato <- function(df) {
    times <- unique(df$time_mandante)
    resultados <- df %>%
        mutate(vitoria = case_when(gols_mandante > gols_visitante ~ time_mandante,
                                   gols_visitante > gols_mandante ~ time_visitante,
                                   gols_mandante == gols_visitante ~ "Empate",
                                   is.na(gols_mandante) ~ NA),
               pontos_mandante = case_when(vitoria == time_mandante ~ 3,
                                           vitoria == "Empate" ~ 1,
                                           vitoria == time_visitante ~ 0,
                                           is.na(vitoria) ~ NA),
               pontos_visitante = case_when(vitoria == time_visitante ~ 3,
                                            vitoria == "Empate" ~ 1,
                                            vitoria == time_mandante ~ 0,
                                            is.na(vitoria) ~ NA))
    pontos_time <- data.frame(times)
    pontos_time$pontos <- 0
    for (i in 1:length(times)) {
        df_mandante <- resultados %>%
            filter(time_mandante == times[i])
        df_visitante <- resultados %>%
            filter(time_visitante == times[i])
        pontos_time$pontos[i] <- sum(c(df_mandante$pontos_mandante, df_visitante$pontos_visitante), na.rm = T)
    }
    pontos_time <- arrange(pontos_time, desc(pontos))
    return(pontos_time)
}

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


simula_jogo <- function(df) {
    gols_m <- rpois(nrow(df %>% filter(rodada >= 26)), (df$theta_m + df$phi_v)/2)
    gols_v <- rpois(nrow(df %>% filter(rodada >= 26)), (df$theta_v + df$phi_m)/2)
    df$gols_mandante[df$rodada >= 26] <- gols_m
    df$gols_visitante[df$rodada >= 26] <- gols_v
    return(df)
}

#====Carregamento dos dados====

dados_bras <- read.csv("~/Downloads/brasileirao_2026.csv")

# Verificando se os resultados condizem com a classificação na 25ª rodada:
dados_bras %>%
    campeonato()

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

taxas_estimadas

df_completo <- cria_df_taxas(dados_bras, taxas_estimadas)

#====Simulações====

set.seed(5050)

t5 <- Sys.time()
simulacao_campeonatos <- replicate(10000, (simula_jogo(df_completo) %>% campeonato)[1, 1])
t6 <- Sys.time()


#====Teste====
library(tidyverse)
library(parallel)

# 1. Separar jogos realizados e futuros
jogos_realizados <- df_completo %>% filter(!is.na(gols_mandante))
jogos_futuros <- df_completo %>% filter(is.na(gols_mandante))

# 2. Calcular a classificação atual UMA ÚNICA VEZ considerando desempates
calcula_tabela <- function(df) {
    df %>%
        mutate(
            pts_m = case_when(gols_mandante > gols_visitante ~ 3, gols_mandante == gols_visitante ~ 1, TRUE ~ 0),
            pts_v = case_when(gols_visitante > gols_mandante ~ 3, gols_mandante == gols_visitante ~ 1, TRUE ~ 0),
            vit_m = ifelse(gols_mandante > gols_visitante, 1, 0),
            vit_v = ifelse(gols_visitante > gols_mandante, 1, 0)
        ) %>%
        pivot_longer(cols = c(time_mandante, time_visitante), names_to = "mando", values_to = "time") %>%
        mutate(
            pontos = ifelse(mando == "time_mandante", pts_m, pts_v),
            vitorias = ifelse(mando == "time_mandante", vit_m, vit_v),
            gols_pro = ifelse(mando == "time_mandante", gols_mandante, gols_visitante),
            gols_contra = ifelse(mando == "time_mandante", gols_visitante, gols_mandante),
            saldo = gols_pro - gols_contra
        ) %>%
        group_by(time) %>%
        summarise(across(c(pontos, vitorias, saldo, gols_pro), sum, na.rm = T), .groups = "drop")
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
t_inicio <- Sys.time()
resultados <- mclapply(1:100000, simula_rodadas_restantes, mc.cores = detectCores() - 1)
t_fim <- Sys.time()
print(t_fim - t_inicio)

lapply(resultados, function(x){x$pontos_todos$pontos[1]}) %>% unlist() %>% hist()
lapply(resultados, function(x){x$pontos_todos$pontos[1]}) %>% unlist() %>% quantile(0.9)

lapply(resultados, function(x){x$pontos_todos$pontos[16]}) %>% unlist() %>% hist()
lapply(resultados, function(x){x$pontos_todos$pontos[16]}) %>% unlist() %>% quantile(0.95)

#====Bootstrap====

B <- 1000 # Número de amostras Bootstrap
boot_theta_m <- matrix(NA, nrow = B, ncol = length(unique(jogos_realizados$time_mandante)),
                       dimnames = list(NULL, unique(jogos_realizados$time_mandante) %>% sort()))
boot_theta_v <- matrix(NA, nrow = B, ncol = length(unique(jogos_realizados$time_mandante)),
                       dimnames = list(NULL, unique(jogos_realizados$time_mandante) %>% sort()))
boot_phi_m <- matrix(NA, nrow = B, ncol = length(unique(jogos_realizados$time_mandante)),
                     dimnames = list(NULL, unique(jogos_realizados$time_mandante) %>% sort()))
boot_phi_v <- matrix(NA, nrow = B, ncol = length(unique(jogos_realizados$time_mandante)),
                     dimnames = list(NULL, unique(jogos_realizados$time_mandante) %>% sort()))

set.seed(5050)

for(i in 1:B) {
    amostra <- jogos_realizados[sample(1:nrow(jogos_realizados), replace = TRUE), ]
    taxas_estimadas_boot <- left_join(
    (amostra %>%
     mutate(time = time_mandante) %>%
     group_by(time) %>%
     summarise(theta_m = mean(gols_mandante, na.rm = T),
               phi_v = mean(gols_visitante, na.rm = T))),
    (amostra %>%
     mutate(time = time_visitante) %>%
     group_by(time) %>%
     summarise(theta_v = mean(gols_visitante, na.rm = T),
               phi_m = mean(gols_mandante, na.rm = T))), by = "time")
    boot_theta_m[i, ] <- taxas_estimadas_boot$theta_m
    boot_theta_v[i, ] <- taxas_estimadas_boot$theta_v
    boot_phi_m[i, ] <- taxas_estimadas_boot$phi_m
    boot_phi_v[i, ] <- taxas_estimadas_boot$phi_v
  # Recalcule as taxas médias de gols para a amostra gerada e armazene
}

ordem_theta_m <- boot_theta_m %>%
    apply(2, median) %>%
    sort() %>%
    names()

stack(as.data.frame(boot_theta_m)) %>%
    mutate(ind = factor(ind, levels = ordem_theta_m)) %>%
    ggplot(aes(x = ind, y = values)) +
    geom_boxplot() +
    geom_point(data = taxas_estimadas %>%
                   mutate(time = factor(time, levels = ordem_theta_v)),
               aes(x = time, y = theta_m), colour = "red") +
    ggpubr::theme_pubclean() +
    labs(x = "Times", y = expression(theta[m])) +
    theme(axis.text.x = element_text(angle = 90))

ordem_theta_v <- boot_theta_v %>%
    apply(2, median) %>%
    sort() %>%
    names()

stack(as.data.frame(boot_theta_v)) %>%
    mutate(ind = factor(ind, levels = ordem_theta_v)) %>%
    ggplot(aes(x = ind, y = values)) +
    geom_boxplot() +
    geom_point(data = taxas_estimadas %>%
                   mutate(time = factor(time, levels = ordem_theta_v)),
               aes(x = time, y = theta_v), colour = "red") +
    ggpubr::theme_pubclean() +
    labs(x = "Times", y = expression(theta[v])) +
    theme(axis.text.x = element_text(angle = 90))

ordem_phi_m <- boot_phi_m %>%
    apply(2, median) %>%
    sort() %>%
    names()

stack(as.data.frame(boot_phi_m)) %>%
    mutate(ind = factor(ind, levels = ordem_phi_m)) %>%
    ggplot(aes(x = ind, y = values)) +
    geom_boxplot() +
    geom_point(data = taxas_estimadas %>%
                   mutate(time = factor(time, levels = ordem_theta_v)),
               aes(x = time, y = phi_m), colour = "red") +
    ggpubr::theme_pubclean() +
    labs(x = "Times", y = expression(Phi[m])) +
    theme(axis.text.x = element_text(angle = 90))

ordem_phi_v <- boot_phi_v %>%
    apply(2, median) %>%
    sort() %>%
    names()

stack(as.data.frame(boot_phi_v)) %>%
    mutate(ind = factor(ind, levels = ordem_phi_v)) %>%
    ggplot(aes(x = ind, y = values)) +
    geom_boxplot() +
    geom_point(data = taxas_estimadas %>%
                   mutate(time = factor(time, levels = ordem_theta_v)),
               aes(x = time, y = phi_v), colour = "red") +
    ggpubr::theme_pubclean() +
    labs(x = "Times", y = expression(Phi[v])) +
    theme(axis.text.x = element_text(angle = 90))

#====Resultados====
left_join(
    data.frame(Time = unique(jogos_futuros$time_mandante) %>% sort()),
    lapply(resultados, function(x){x$campeao}) %>% unlist() %>% table() %>% as.data.frame() %>% mutate(Freq = Freq/100000) %>% select("Time" = ".","Freq")) %>%
    arrange(desc(Freq)) %>%
    mutate(Freq = case_when(is.na(Freq) ~ 0,
                            .default = Freq),
           Freq = scales::percent(Freq)) %>%
    select(Time, "Prob_Camp" = Freq)

left_join(
    data.frame(Time = unique(jogos_futuros$time_mandante) %>% sort()),
    lapply(resultados, function(x){x$rebaixados}) %>% unlist() %>% table() %>% as.data.frame() %>% mutate(Freq = Freq/400000) %>% select("Time" = ".","Freq")) %>%
    arrange(desc(Freq)) %>%
    mutate(Freq = case_when(is.na(Freq) ~ 0,
                            .default = Freq),
           Freq = scales::percent(Freq)) %>%
    select(Time, "Prob_Reb" = Freq)

#==== Pacotes necessários ====
library(tidyverse)
library(parallel)
library(purrr)

#==== 1. Leitura e Separação dos Dados ====
# Substitua o caminho pelo local correto do seu arquivo
dados_bras <- read.csv("~/Downloads/brasileirao_2026.csv")

jogos_realizados <- dados_bras %>% filter(!is.na(gols_mandante))
jogos_futuros <- dados_bras %>% filter(is.na(gols_mandante))

#==== 2. Cálculo das Taxas Variáveis via EWMA ====
# O parâmetro alpha controla o "esquecimento" do modelo.
# alpha = 0.2 significa que o jogo atual tem 20% de peso na nova taxa, e o histórico 80%.
alpha <- 0.20

# Função auxiliar para calcular o EWMA iterativamente usando purrr::accumulate
calcula_ewma <- function(vetor_gols) {
  accumulate(vetor_gols, ~ alpha * .y + (1 - alpha) * .x, .init = vetor_gols[1])[-1]
}

# Criamos o histórico linear, calculamos o EWMA rodada a rodada e extraímos a última taxa de cada time
taxas_finais_ewma <- jogos_realizados %>%
  pivot_longer(
    cols = c(time_mandante, time_visitante),
    names_to = "mando",
    values_to = "time"
  ) %>%
  mutate(
    gols_feitos = ifelse(mando == "time_mandante", gols_mandante, gols_visitante),
    gols_sofridos = ifelse(mando == "time_mandante", gols_visitante, gols_mandante)
  ) %>%
  arrange(time, rodada) %>%
  group_by(time) %>%
  mutate(
    theta_hist = calcula_ewma(gols_feitos),
    phi_hist = calcula_ewma(gols_sofridos)
  ) %>%
  # Pegamos apenas o valor da taxa na última rodada jogada pelo time para projetar o futuro
  slice_tail(n = 1) %>%
  select(time, theta_ewma = theta_hist, phi_ewma = phi_hist) %>%
  ungroup()

#==== 3. Preparação da Tabela Base e Jogos Futuros ====
# Cruzamos as taxas atualizadas com a tabela de jogos ainda não realizados
jogos_futuros_ewma <- jogos_futuros %>%
    select(-c(theta_m, theta_v, phi_m, phi_v)) %>%
  left_join(
    taxas_finais_ewma %>% rename(time_mandante = time, theta_m = theta_ewma, phi_m = phi_ewma),
    by = "time_mandante"
  ) %>%
  left_join(
    taxas_finais_ewma %>% rename(time_visitante = time, theta_v = theta_ewma, phi_v = phi_ewma),
    by = "time_visitante"
  )

calcula_tabela <- function(df) {
  df %>%
    mutate(
      pts_m = case_when(gols_mandante > gols_visitante ~ 3, gols_mandante == gols_visitante ~ 1, TRUE ~ 0),
      pts_v = case_when(gols_visitante > gols_mandante ~ 3, gols_mandante == gols_visitante ~ 1, TRUE ~ 0),
      vit_m = ifelse(gols_mandante > gols_visitante, 1, 0),
      vit_v = ifelse(gols_visitante > gols_mandante, 1, 0)
    ) %>%
    pivot_longer(cols = c(time_mandante, time_visitante), names_to = "mando", values_to = "time") %>%
    mutate(
      pontos = ifelse(mando == "time_mandante", pts_m, pts_v),
      vitorias = ifelse(mando == "time_mandante", vit_m, vit_v),
      gols_pro = ifelse(mando == "time_mandante", gols_mandante, gols_visitante),
      gols_contra = ifelse(mando == "time_mandante", gols_visitante, gols_mandante),
      saldo = gols_pro - gols_contra
    ) %>%
    group_by(time) %>%
    summarise(across(c(pontos, vitorias, saldo, gols_pro), sum), .groups = "drop")
}

# Computa a classificação real do campeonato até o momento exato da pausa
tabela_base <- calcula_tabela(jogos_realizados)

#==== 4. Função de Simulação ====
simula_campeonato_ewma <- function(...) {
  sim_futuros <- jogos_futuros_ewma %>%
    mutate(
        gols_mandante = rpois(n(), (theta_m + phi_v)/2),
        gols_visitante = rpois(n(), (theta_v + phi_m)/2)
    )
  tabela_simulada <- calcula_tabela(sim_futuros)
  # Junta o passado real com o futuro simulado e aplica os critérios de desempate
  tabela_final <- bind_rows(tabela_base, tabela_simulada) %>%
      group_by(time) %>%
      summarise(across(everything(), sum), .groups = "drop") %>%
      arrange(desc(pontos), desc(vitorias), desc(saldo), desc(gols_pro)) %>%
      mutate(posicao = row_number())
  list(
      campeao = tabela_final$time[1],
      pontos_campeao = tabela_final$pontos[1],
      vice = tabela_final$time[2],
      pontos_vice = tabela_final$pontos[2],
      rebaixados = tabela_final$time[17:20],
      pontos_todos = tabela_final %>% select(time, pontos, posicao)
  )
}

#==== 5. Execução por Monte Carlo ====
set.seed(5050)

# Aproveita os múltiplos núcleos de processamento para rodar as 10.000 iterações
n_simulacoes <- 100000
t_inicio <- Sys.time()
resultados_mc <- mclapply(1:n_simulacoes, simula_campeonato_ewma, mc.cores = detectCores() - 2)
t_fim <- Sys.time()

cat("Tempo de execução das simulações:", round(difftime(t_fim, t_inicio, units = "secs"), 2), "segundos\n")

left_join(
    data.frame(Time = unique(jogos_futuros$time_mandante) %>% sort()),
    lapply(resultados_mc, function(x){x$campeao}) %>% unlist() %>% table() %>% as.data.frame() %>% mutate(Freq = Freq/100000) %>% select("Time" = ".","Freq")) %>%
    arrange(desc(Freq)) %>%
    mutate(Freq = case_when(is.na(Freq) ~ 0,
                            .default = Freq),
           Freq = scales::percent(Freq)) %>%
    select(Time, "Prob_Camp" = Freq)

left_join(
    data.frame(Time = unique(jogos_futuros$time_mandante) %>% sort()),
    lapply(resultados_mc, function(x){x$rebaixados}) %>% unlist() %>% table() %>% as.data.frame() %>% mutate(Freq = Freq/400000) %>% select("Time" = ".","Freq")) %>%
    arrange(desc(Freq)) %>%
    mutate(Freq = case_when(is.na(Freq) ~ 0,
                            .default = Freq),
           Freq = scales::percent(Freq)) %>%
    select(Time, "Prob_Reb" = Freq)

#====Gráficos====

# Evolução do theta via EWMA:
jogos_realizados %>%
    pivot_longer(
        cols = c(time_mandante, time_visitante),
        names_to = "mando",
        values_to = "time"
    ) %>%
    mutate(
        gols_feitos = ifelse(mando == "time_mandante", gols_mandante, gols_visitante),
        gols_sofridos = ifelse(mando == "time_mandante", gols_visitante, gols_mandante)
    ) %>%
    arrange(time, rodada) %>%
    group_by(time) %>%
    mutate(
        theta_hist = calcula_ewma(gols_feitos),
        phi_hist = calcula_ewma(gols_sofridos)
    ) %>%
    filter(time %in% c("Athletico-PR", "Coritiba")) %>%
    ggplot(aes(x = rodada, y = theta_hist, colour = time)) +
    geom_point() +
    geom_line() +
    theme_bw()

set.seed(123)
jogos_realizados %>%
    pivot_longer(
        cols = c(time_mandante, time_visitante),
        names_to = "mando",
        values_to = "time"
    ) %>%
    mutate(
        gols_feitos = ifelse(mando == "time_mandante", gols_mandante, gols_visitante),
        gols_sofridos = ifelse(mando == "time_mandante", gols_visitante, gols_mandante)
    ) %>%
    arrange(time, rodada) %>%
bind_rows(
(jogos_futuros_ewma %>%
 mutate(
     gols_mandante = rpois(n(), (theta_m + phi_v)/2),
     gols_visitante = rpois(n(), (theta_v + phi_m)/2)
 ) %>% pivot_longer(
                cols = c(time_mandante, time_visitante),
                names_to = "mando",
                values_to = "time"
            ) %>%
  mutate(
    gols_feitos = ifelse(mando == "time_mandante", gols_mandante, gols_visitante),
    gols_sofridos = ifelse(mando == "time_mandante", gols_visitante, gols_mandante)
  ) %>%
  arrange(time, rodada))) %>%
  group_by(time) %>%
  mutate(
    theta_hist = calcula_ewma(gols_feitos),
    phi_hist = calcula_ewma(gols_sofridos)
  ) %>%
    filter(time %in% c("Athletico-PR", "Coritiba")) %>%
    mutate(rod = case_when(rodada >= 26 ~ "Simulado",
                           .default = "Real")) %>%
    ggplot(aes(x = rodada, y = theta_hist, colour = time, shape = rod)) +
    geom_point() +
    geom_line() +
    theme_bw()
