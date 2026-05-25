#----Pacotes----

library(bayesianETAS)

library(readxl)

library(ETAS)

library(tidyverse)

library(INLA)

library(inlabru)

source("Hawkes_process_tutorials/code_for_paper/hawkes_functions.R")

source("Hawkes_process_tutorials/code_for_paper/hawkes_utilities.R")

#----Base de dados----

# Transformações necessárias nos dados:
peru <- read_xlsx("Dados/Inst_Peru.xlsx") %>%
    select(date = `fecha UTC`, time = `hora UTC`,
           lat = `latitud (º)`, long = `longitud (º)`,
           prof = `profundidad (km)`,
           mag = `magnitud (M)`) %>%
    mutate(
        date = as.Date(date),
        lat = as.numeric(lat),
        long = as.numeric(long),
        prof = as.numeric(prof),
        mag = as.numeric(mag)
    ) %>%
    filter(mag >= 4)

# Catálogo para pacote ETAS:
peru.quakes <- catalog(peru)

# Gráfico padrão do pacote ETAS:
plot(peru.quakes)

# Animação com os últimos 1.000 terremotos:
stpp::animation(stpp::as.3dpoints(tail(peru, 1000)$lat,
                                  tail(peru, 1000)$long,
                                  tail(peru, 1000)$date),
                runtime = 10, cex = tail(peru, 1000)$mag/3)

stpp::animation(stpp::as.3dpoints(dados_simulados$x,
                                  dados_simulados$y,
                                  dados_simulados$t),
                runtime = 10, cex = dados_simulados$m/3)

# Not run:
aj_etas <- etas(peru.quakes)

aj_bayes <- maxLikelihoodETAS(as.numeric(ymd_hms(str_c(peru$date, "-", peru$time))), peru$mag, 4,
                              max(as.numeric(ymd_hms(str_c(peru$date, "-", peru$time)))))

#----INLA----

# Magnitude mínima:
M0 = 3.99

# Tempos mínimos e máximos:
T1 = 0
T2 = 24113

# Importar os dados para formato especificado:
peru_inla <- read_xlsx(path = "Dados/Inst_Peru.xlsx") %>%
    select(date = `fecha UTC`, time = `hora UTC`,
           lat = `latitud (º)`, long = `longitud (º)`,
           prof = `profundidad (km)`,
           mag = `magnitud (M)`) %>%
    mutate(
        time_date = as.POSIXct(paste0(date, time)),
        time.diff = as.numeric(difftime(time_date, min(time_date) - 1, units = 'days')),
        lat = as.numeric(lat),
        long = as.numeric(long),
        prof = as.numeric(prof),
        mag = as.numeric(mag),
        mag.class = cut(mag, breaks = c(M0, 5, 7))
    ) %>%
    arrange(time_date)

# Dados para Inlabru:
data.bru <- data.frame(ts = peru_inla$time.diff,
                       magnitudes = peru_inla$mag) %>%
  mutate(idx.p = 1:nrow(peru_inla))

# Valores iniciais para theta:
th.init <- list(th.mu = 0.5,
                th.K = 0.5,
                th.alpha = -2,
                th.c = -2,
                th.p = -2)

# Funções de ligação para prioris:
link.f.be <- list(mu = \(x) gamma.t(x, 0.1, 0.1),
                  K = \(x) loggaus.t(x, -1, 2),
                  alpha = \(x) unif.t(x, 0, 10),
                  c_ = \(x) unif.t(x, 0, 10),
                  p = \(x) unif.t(x, 1, 10))

# Opções para Inlabru:
bru.opt.list <- list(bru_verbose = 4, # Tipo de output visual
                     bru_max_iter = 100, # Número máximo de iterações
                     #bru_method = list(max_step = 0.5, rel_tol = 0.01),
                     inla.mode = 'experimental', # Tipo de algoritmo usado
                     bru_initial = th.init) # Parâmetros iniciais

# Ajuste usando Inlabru:
aj_inla <- Hawkes.bru(sample.s = data.bru, # Dados
                      M0 = M0, # Magnitude mínima
                      T1 = 0, T2 = T2, # Domínio temporal
                      link.functions = link.f.be, # Funções de ligação
                      coef.t. = 1, # binning parameter (delta)
                      delta.t. = 0.1, # binning parameter (Delta)
                      N.max. = 3, # binning parameter (n.max)
                      bru.opt = bru.opt.list)

# Predição INLA:
predict(aj_inla,
        data.frame(),
        ~ lambda.N(th.mu, th.K, th.alpha, th.c, th.p,
                   T1, T2, M0,
                   data.bru,
                   link.f.be))

predict(aj_inla, data.frame(),
                         ~ data.frame(N = 24000:26000,
                                      pdf = dpois(24000:26000,
                                                  lambda.N(th.mu, th.K,
                                                           th.alpha, th.c, th.p,
                                                           T1, T2, M0,
                                                           data.bru,
                                                           link.f.be)))) %>%
    ggplot(aes(x = N, y = mean)) +
    geom_line(color = 'darkblue') +
    geom_ribbon(aes(xmax = N, xmin = N, ymin = q0.025, ymax = q0.975), alpha = 0.2,
                fill = 'blue') +
    geom_vline(xintercept = nrow(data.bru), linetype = 3) +
    geom_line(data = data.frame(x = 24000:26000,
                                y = dpois(24000:26000, 967.6161)),
              aes(x,y), color = 'red', linetype = 2) +
    ylab('pdf')

#----Algoritmos de Simulação----

# Simulação de Processo de Hawkes Espaço-Temporal Marcado (Bernabeu et. al. 2024):

simulate_hawkes_marked <- function(
                                   T_max = 100,           # Tempo máximo
                                   X_lim = c(0, 10),      # Limites em X
                                   Y_lim = c(0, 10),      # Limites em Y
                                   mu = 0.5,              # Intensidade de base
                                   # Parâmetros do Kernel g:
                                   k0 = 0.4,              # Produtividade base (o quanto um terremoto afeta a criação de outros)
                                   alpha = 0.5,           # Influência da magnitude na produtividade
                                   beta_time = 1.0,       # Decaimento temporal (kernel exponencial)
                                   sigma_space = 0.5,     # Desvio padrão espacial (kernel gaussiano)
                                   b_mag = 1.0,           # Parâmetro da distribuição exponencial das magnitudes
                                   m_min = 2.0            # Magnitude mínima a ser considerada
                                   ) {
    # Lista para armazenar todos os eventos:
    all_events <- data.frame(t=numeric(), x=numeric(), y=numeric(), m=numeric(), gen=numeric())
    #----Geração dos Eventos Pais----
    # Área da região:
    area <- (X_lim[2] - X_lim[1]) * (Y_lim[2] - Y_lim[1])
    # Número esperado de pais (Poisson):
    n_parents <- rpois(1, lambda = mu * area * T_max)
    if(n_parents > 0) {
        # Criação dos eventos pais:
        parents <- data.frame(
            t = runif(n_parents, 0, T_max),            # Tempos distribuídos uniformemente
            x = runif(n_parents, X_lim[1], X_lim[2]),  # Posições distribuídas uniformemente
            y = runif(n_parents, Y_lim[1], Y_lim[2]),  # Posições distribuídas uniformemente
            m = rexp(n_parents, rate = b_mag) + m_min, # Marcas (Magnitude)
            gen = 0                                    # Geração 0
        )
        # Adição à lista com todos os eventos:
        all_events <- rbind(all_events, parents)
        # Controle das gerações:
        current_gen_events <- parents
    } else {
        # Controle para caso nenhum evento seja gerado
        return(NULL)
    }
    #-----Geração de Descendentes----
    # Inicializa o processo na geração 0:
    generation <- 0
    # Loop enquanto houver eventos na geração atual capazes de gerar filhos:
    while(nrow(current_gen_events) > 0) {
        # Atualiza a geração:
        generation <- generation + 1
        # Cria o data.frame com os próximos eventos:
        next_gen_events <- data.frame()
        # Itera sobre cada evento da geração atual:
        for(i in 1:nrow(current_gen_events)) {
            # Adiciona os pais:
            parent <- current_gen_events[i,]
            # Calcula a produtividade baseada na magnitude:
            # N ~ Poisson(k0 * exp(alpha * (m - m_min)))
            # Número esperado de descendentes:
            expected_offspring <- k0 * exp(alpha * (parent$m - m_min))
            # Gera o número de descendentes:
            n_offspring <- rpois(1, expected_offspring)
            if(n_offspring > 0) {
                # Gera tempos dos descendentes, seguindo distribuição exponencial:
                dt <- rexp(n_offspring, rate = beta_time)
                t_new <- parent$t + dt
                # Filtra eventos que passam do T_max:
                valid_idx <- t_new <= T_max
                n_valid <- sum(valid_idx)
                if(n_valid > 0) {
                    t_new <- t_new[valid_idx]
                    # Gera distribuição isotrópica no espaço:
                    x_new <- rnorm(n_valid, mean = parent$x, sd = sigma_space)
                    y_new <- rnorm(n_valid, mean = parent$y, sd = sigma_space)
                    # Gera marcas/magnitudes para os filhos:
                    m_new <- rexp(n_valid, rate = b_mag) + m_min
                    # Cria dataframe temporário:
                    offspring <- data.frame(
                        t = t_new,
                        x = x_new,
                        y = y_new,
                        m = m_new,
                        gen = generation
                    )
                    # Manter apenas eventos dentro da janela espacial:
                    in_window <- offspring$x >= X_lim[1] & offspring$x <= X_lim[2] &
                        offspring$y >= Y_lim[1] & offspring$y <= Y_lim[2]
                    if(sum(in_window) > 0) {
                        # Adiciona os eventos desta geração com os das gerações anteriores:
                        next_gen_events <- rbind(next_gen_events, offspring[in_window, ])
                    }
                }
            }
        }
        # Atualiza lista principal e preparar para próxima iteração:
        if(nrow(next_gen_events) > 0) {
            all_events <- rbind(all_events, next_gen_events)
            current_gen_events <- next_gen_events
        } else {
            # Encerra o loop, quando não há mais eventos criados na geração:
            current_gen_events <- data.frame() # Encerra o loop
        }
        # Segurança para evitar loops infinitos em casos supercríticos:
        if(generation > 100){
            break
        }
    }
    # Ordena os eventos por tempo:
    all_events <- all_events[order(all_events$t), ]
    # Retorna o data.frame com todos os eventos ordenados:
    return(all_events)
}

# Simulação de Processo de Hawkes Espaço-Temporal Baseado em ETAS (Ogata 1998):

simulate_hawkes_etas <- function(
                                 T_max = 100,           # Tempo máximo:
                                 X_lim = c(-10, 10),    # Limites em X
                                 Y_lim = c(-10, 10),    # Limites em Y
                                 mu = 0.05,             # Intensidade de base
                                 # Parâmetros do Kernel g:
                                 k0 = 0.1,              # Constante de escala
                                 c_par = 0.5,           # Parâmetro de deslocamento temporal
                                 p_par = 2.5,           # Expoente de decaimento temporal (p > 1 para convergir)
                                 d_par = 0.2,           # Escala espacial base
                                 alpha = 1.0,           # Sensibilidade à magnitude
                                 m_min = 2.0,           # Magnitude mínima a ser considerada
                                 b_mag = 1.5            # Parâmetro da distribuição exponencial das magnitudes
                                 ) {
    #----Funções Auxiliares----
    # Integral Temporal (0 a Infinito) de (t+c)^-p
    get_time_integral <- function(c, p) {
        return (c^(1-p) / (p-1))
    }
    # Integral Espacial de exp(-0.5 * r^2 / D)
    get_space_integral <- function(D) {
        return (2 * pi * D)
    }
    # Amostragem Temporal (Lei de Omori para Pareto)
    r_omori <- function(n, c, p) {
        u <- runif(n)
        return (c * ((1 - u)^(1 / (1 - p)) - 1))
    }
    # Lista para armazenar todos os eventos:
    all_events <- data.frame(t=numeric(), x=numeric(), y=numeric(), m=numeric(), gen=numeric())
    #----Geração dos eventos pais---
    # Área da região:
    area <- (X_lim[2] - X_lim[1]) * (Y_lim[2] - Y_lim[1])
    # Número de eventos pais gerados:
    n_parents <- rpois(1, lambda = mu * area * T_max)
    if(n_parents > 0) {
        # Criação dos eventos pais:
        current_gen <- data.frame(
            t = runif(n_parents, 0, T_max),              # Tempos uniformemente distribuídos
            x = runif(n_parents, X_lim[1], X_lim[2]),    # Localizações uniformemente distribuídas
            y = runif(n_parents, Y_lim[1], Y_lim[2]),    # Localizações uniformemente distribuídas
            m = rexp(n_parents, rate = b_mag) + m_min,   # Magnitudes exponenciais
            gen = 0
        )
        all_events <- rbind(all_events, current_gen)
    } else {
        # Caso onde nenhum evento tenha sido gerado:
        return(NULL)
    }
    #----Geração de Descendentes----
    # Inicializa a contagem de gerações:
    generation <- 0
    while(nrow(current_gen) > 0) {
        generation <- generation + 1
        next_gen <- data.frame()
        for(i in 1:nrow(current_gen)) {
            # Adiciona os eventos pais:
            parent <- current_gen[i,]
            # Calcula o termo D dependente da magnitude:
            D_val <- d_par * exp(alpha * (parent$m - m_min))
            # Calcula a produtividade (Número esperado de filhos):
            # N_esp = k0 * Int_Tempo * Int_Espaço
            int_t <- get_time_integral(c_par, p_par)
            int_s <- get_space_integral(D_val)
            expected_offspring <- k0 * int_t * int_s
            n_offspring <- rpois(1, expected_offspring)
            if(n_offspring > 0) {
                # Gera tempos aleatórios (seguindo Lei de Omori):
                dt <- r_omori(n_offspring, c_par, p_par)
                t_new <- parent$t + dt
                # Filtro temporal:
                valid <- t_new <= T_max
                if(sum(valid) > 0) {
                    n_val <- sum(valid)
                    t_new <- t_new[valid]
                    # Gera distribuição isométrica no espaço (Gaussiana com variância D_val):
                    sigma <- sqrt(D_val)
                    x_new <- rnorm(n_val, mean = parent$x, sd = sigma)
                    y_new <- rnorm(n_val, mean = parent$y, sd = sigma)
                    # Gera magnitudes dos descendentes:
                    m_new <- rexp(n_val, rate = b_mag) + m_min
                    offspring <- data.frame(t = t_new, x = x_new, y = y_new, m = m_new, gen = generation)
                    # Filtro espacial:
                    in_window <- offspring$x >= X_lim[1] & offspring$x <= X_lim[2] &
                        offspring$y >= Y_lim[1] & offspring$y <= Y_lim[2]
                    if(sum(in_window) > 0) {
                        # Adiciona os descendentes à lista de eventos:
                        next_gen <- rbind(next_gen, offspring[in_window,])
                    }
                }
            }
        }
        if(nrow(next_gen) > 0) {
            # Cria lista com todos os eventos:
            all_events <- rbind(all_events, next_gen)
            current_gen <- next_gen
        } else {
            break
        }
        # Segurança para evitar loops infinitos em casos supercríticos:
        if(generation > 100){
            break
        }
    }
    # Retorna todos os eventos, ordenados pelo tempo:
    return(all_events[order(all_events$t), ])
}

#----Exemplos----

# Primeiro simulador:
set.seed(5050)

dados_simulados <- simulate_hawkes_marked(
  T_max = 12,            # Tempo máximo de 12 anos
  mu = 0.1,              # Intensidade base
  k0 = 0.5,              # Produtividade base
  b_mag = 1.5,           # Beta da distribuição exponencial das magnitudes
  alpha = 0.7,           # Modificador de intensidade
  sigma_space = 0.3,     # Variabilidade espacial
  m_min = 3              # Magnitude mínima
)

ggplot(dados_simulados, aes(x = x, y = y)) +
    geom_point(aes(color = t, size = m), alpha = 0.7) +
    scale_color_viridis_c(option = "plasma", name = "Tempo") +
    scale_size_continuous(name = "Magnitude") +
    theme_minimal() +
    labs(title = "Simulação Hawkes Espaço-Temporal Marcado",
         subtitle = "Kernel baseado em Bernabeu et. al. 2024") +
    coord_fixed()

ggplot(dados_simulados, aes(x = x, y = y)) +
    geom_point(aes(color = gen, size = m), alpha = 0.7) +
    scale_color_viridis_c(option = "plasma", name = "Geração") +
    scale_size_continuous(name = "Magnitude") +
    theme_minimal() +
    labs(title = "Simulação Hawkes Espaço-Temporal Marcado",
         subtitle = "Kernel baseado em Bernabeu et. al. 2024") +
    coord_fixed()

# Segundo simulador:
dados_simulados2 <- simulate_hawkes_etas(
    T_max = 12,         # Tempo máximo de 12 anos
    X_lim = c(0, 10),   # Limite da caixa em X
    Y_lim = c(0, 10),   # Limite da caixa em Y
    mu = 0.1,           # Intensidade base
    k0 = 0.5,           # Produtividade base
    alpha = 0.7,        # Modificador de intensidade
    b_mag = 1.5,        # Beta da distribuição exponencial das magnitudes
    p_par = 1.1,        # Decaimento temporal
    m_min = 3,          # Magnitude mínima
    d_par = 0.1         # Escala espacial
)

ggplot(dados_simulados2, aes(x = x, y = y)) +
    geom_point(aes(color = t, size = m), alpha = 0.7) +
    scale_color_viridis_c(option = "plasma", name = "Tempo") +
    scale_size_continuous(name = "Magnitude") +
    theme_minimal() +
    labs(title = "Simulação Hawkes Espaço-Temporal Marcado",
         subtitle = "Kernel baseado em Ogata 1998") +
    coord_fixed()

ggplot(dados_simulados2, aes(x = x, y = y)) +
    geom_point(aes(color = gen, size = m), alpha = 0.7) +
    scale_color_viridis_c(option = "plasma", name = "Geração") +
    scale_size_continuous(name = "Magnitude") +
    theme_minimal() +
    labs(title = "Simulação Hawkes Espaço-Temporal Marcado",
         subtitle = "Kernel baseado em Ogata 1998") +
    coord_fixed()

#----Comparação----

# Quantidade de eventos por simulação:
cat("Total de eventos em 1:", nrow(dados_simulados), "\n",
    "Total de eventos em 2:", nrow(dados_simulados2), "\n")

# Contagem de eventos por geração:
table(dados_simulados$gen)
table(dados_simulados2$gen)

#----Animações----

anim <- ggplot(dados_simulados, aes(x = x, y = y)) +
    geom_point(aes(color = gen, size = m, group = t %/% (1/12)), alpha = 0.7) +
    scale_color_viridis_c(option = "plasma", name = "Geração") +
    scale_size_continuous(name = "Magnitude") +
    theme_minimal() +
    coord_fixed() +
    gganimate::transition_states(t %/% (1/12)) +
    ggtitle('Mês : {closest_state}')

gganimate::animate(anim)

anim2 <- ggplot(dados_simulados2, aes(x = x, y = y)) +
    geom_point(aes(color = gen, size = m, group = t %/% (1/4)), alpha = 0.7) +
    scale_color_viridis_c(option = "plasma", name = "Geração") +
    scale_size_continuous(name = "Magnitude") +
    theme_minimal() +
    coord_fixed() +
    gganimate::transition_states(t %/% (1/4)) +
    ggtitle('Trimestre : {closest_state}')

gganimate::animate(anim2)


#----Ajustes nas Simulações----

simulate_hawkes_etas_jss <- function(
                                     # Parâmetros gerais:
                                     T_max = 100,           # Tempo máximo
                                     X_lim = c(-10, 10),    # Limites em X
                                     Y_lim = c(-10, 10),    # Limites em Y
                                     mu = 0.5,              # Taxa de background base
                                     u_fn = NULL,           # Função não-homogênea u(x,y) que retorna a intensidade espacial
                                     u_max = 1.0,           # Valor máximo da função u_fn na janela para aceitação/rejeição
                                     # Parâmetros ETAS:
                                     A = 0.1,               # Produtividade
                                     alpha = 1.05,          # Sensibilidade da produtividade à magnitude
                                     c_par = 0.05,          # Deslocamento temporal
                                     p_par = 1.2,           # Decaimento temporal (p > 1)
                                     D_par = 0.02,          # Escala espacial base
                                     gamma = 0.5,           # Efeito da magnitude na escala espacial
                                     q_par = 1.5,           # Decaimento espacial (q > 1)
                                     beta = 1.5,            # Parâmetro da exponencial para distribuição das magnitudes
                                     m_min = 2.0            # Magnitude mínima (m0)
                                     ){
    # Lista para armazenar todos os eventos:
    all_events <- data.frame(t = numeric(), x = numeric(), y = numeric(), m = numeric(), gen = numeric())
    area <- (X_lim[2] - X_lim[1]) * (Y_lim[2] - Y_lim[1])
    # ---- 1. Geração dos Eventos Pais (Background Não-Homogêneo) ----
    # Algoritmo de Thinning (Aceitação-Rejeição)
    # Caso seja especificada função de intensidade:
    if (!is.null(u_fn)) {
        lambda_max <- mu * u_max * area * T_max
        n_candidates <- rpois(1, lambda = lambda_max)
        if (n_candidates > 0) {
            c_x <- runif(n_candidates, X_lim[1], X_lim[2])
            c_y <- runif(n_candidates, Y_lim[1], Y_lim[2])
            # Probabilidade de aceitar o evento baseada em u(x,y)
            prob_accept <- u_fn(c_x, c_y) / u_max
            accepted <- runif(n_candidates) < prob_accept
            n_parents <- sum(accepted)
            if (n_parents > 0) {
                parents <- data.frame(
                    t = runif(n_parents, 0, T_max),
                    x = c_x[accepted],
                    y = c_y[accepted],
                    m = rexp(n_parents, rate = beta) + m_min,
                    gen = 0
                )
                all_events <- rbind(all_events, parents)
                current_gen <- parents
            } else {
                return(NULL)
            }
        } else {
            return(NULL)
        }
    } else {
        # Processo de Poisson Homogêneo se u_fn não for fornecido
        n_parents <- rpois(1, lambda = mu * area * T_max)
        if (n_parents > 0) {
            parents <- data.frame(
                t = runif(n_parents, 0, T_max),
                x = runif(n_parents, X_lim[1], X_lim[2]),
                y = runif(n_parents, Y_lim[1], Y_lim[2]),
                m = rexp(n_parents, rate = beta) + m_min,
                gen = 0
            )
            all_events <- rbind(all_events, parents)
            current_gen <- parents
        } else {
            return(NULL)
        }
    }
    # ---- 2. Geração de Descendentes ----
    generation <- 0
    while(nrow(current_gen) > 0) {
        generation <- generation + 1
        next_gen <- data.frame()
        for(i in 1:nrow(current_gen)) {
            parent <- current_gen[i,]
            # Produtividade: N ~ Poisson( A * exp(alpha * (m - m_min)) )
            expected_offspring <- A * exp(alpha * (parent$m - m_min))
            n_offspring <- rpois(1, expected_offspring)
            if(n_offspring > 0) {
                # 2a. Geração de Tempos (Inversa da CDF de Omori modificada)
                U_t <- runif(n_offspring)
                dt <- c_par * ((1 - U_t)^(-1 / (p_par - 1)) - 1)
                t_new <- parent$t + dt
                valid_t <- t_new <= T_max
                n_val <- sum(valid_t)
                if(n_val > 0) {
                    t_new <- t_new[valid_t]
                    # 2b. Geração do Espaço (Inversa da CDF da distribuição pareto)
                    # Variância espacial dependente da magnitude
                    S_val <- D_par * exp(gamma * (parent$m - m_min))
                    U_s <- runif(n_val)
                    R_squared <- S_val * ((1 - U_s)^(-1 / (q_par - 1)) - 1)
                    R <- sqrt(R_squared)
                    theta <- runif(n_val, 0, 2 * pi)
                    x_new <- parent$x + R * cos(theta)
                    y_new <- parent$y + R * sin(theta)
                    # 2c. Marcas / Magnitudes
                    m_new <- rexp(n_val, rate = beta) + m_min
                    offspring <- data.frame(
                        t = t_new,
                        x = x_new,
                        y = y_new,
                        m = m_new,
                        gen = generation
                    )
                    # Filtro espacial para manter na janela de observação
                    in_window <- offspring$x >= X_lim[1] & offspring$x <= X_lim[2] &
                        offspring$y >= Y_lim[1] & offspring$y <= Y_lim[2]
                    if(sum(in_window) > 0) {
                        next_gen <- rbind(next_gen, offspring[in_window,])
                    }
                }
            }
        }
        if(nrow(next_gen) > 0) {
            all_events <- rbind(all_events, next_gen)
            current_gen <- next_gen
        } else {
            break
        }
        # Prevenção de loops infinitos para regimes explosivos/supercríticos
        if(generation > 100){
            warning("Regime supercrítico detectado (mais de 100 gerações). Interrompendo.")
            break
        }
    }
    return(all_events[order(all_events$t), ])
}

# Definindo uma função de background não-homogêneo:
my_u_fn <- function(x, y) {
    return(1)
  # return(exp(-0.05 * (x^2 + y^2)))
}

val_x <- seq(-10, 10, length.out = 100)
val_y <- seq(-10, 10, length.out = 100)

grid <- expand.grid(val_x, val_y)

val_z <- apply(grid, 1, function(x){my_u_fn(x[1], x[2])})

fields::image.plot(val_x,
                   val_y,
                   matrix(val_z,
                          nrow = 100))

set.seed(5050)

sim_dados_homog <- simulate_hawkes_etas_jss(
  T_max = 20,
  X_lim = c(-10, 10), Y_lim = c(-10, 10),
  mu = 1.0, u_fn = my_u_fn, u_max = 1.0,
  A = 0.05, alpha = 1, c_par = 0.1, p_par = 1.3,
  D_par = 0.05, gamma = 1, q_par = 2, beta = 2, m_min = 3
)

catalogo_teste <- catalog(
    sim_dados_homog %>%
    mutate(t = dyears(t),
           data = POSIXct(1) + t,
           long = x, lat = y,
           mag = m,
           date = date(data),
           time = format(data, format = "%H:%M:%S"))
)

plot(catalogo_teste)

# Não roda:
aj_teste_homog_1 <- etas(catalogo_teste, nthreads = 8)

aj_teste_homog_2 <- etas(catalogo_teste,
                         param0 = c(mu = 1, A = 0.02,
                                    ),
                         nthreads = 8)


plot(aj_teste)

rates(aj_teste)
