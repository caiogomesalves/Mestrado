
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

simulate_hawkes_etas_jss <- function(
  T_max = 100,           # Tempo máximo
  X_lim = c(-10, 10),    # Limites em X
  Y_lim = c(-10, 10),    # Limites em Y
  mu = 0.5,              # Taxa de background base
  u_fn = NULL,           # Função não-homogênea u(x,y) que retorna a intensidade espacial
  u_max = 1.0,           # Valor máximo da função u_fn na janela para aceitação/rejeição
  # Parâmetros ETAS (JSS):
  A = 0.1,               # Produtividade
  alpha = 1.05,          # Sensibilidade da produtividade à magnitude
  c_par = 0.05,          # Deslocamento temporal
  p_par = 1.2,           # Decaimento temporal (p > 1)
  D_par = 0.02,          # Escala espacial base
  gamma = 0.5,           # Efeito da magnitude na escala espacial
  q_par = 1.5,           # Decaimento espacial (q > 1)
  beta = 1.5,            # Parâmetro da exponencial para as magnitudes
  m_min = 2.0            # Magnitude mínima (m0)
) {
  # Lista para armazenar todos os eventos:
  all_events <- data.frame(t=numeric(), x=numeric(), y=numeric(), m=numeric(), gen=numeric())
  area <- (X_lim[2] - X_lim[1]) * (Y_lim[2] - Y_lim[1])
  # ---- 1. Geração dos Eventos Pais (Background Não-Homogêneo) ----
  # Usando Thinning (Aceitação-Rejeição)
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
          # 2b. Geração do Espaço (Inversa da CDF da distribuição pareto/bivariada)
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

# --- Exemplo de Uso ---
# Definindo uma função de background não-homogêneo (ex: maior intensidade no centro)
my_u_fn <- function(x, y) {
  return(exp(-0.1 * (x^2 + y^2)))
}

# Simulando
sim_data <- simulate_hawkes_etas_jss(
  T_max = 500,
  X_lim = c(-10, 10), Y_lim = c(-10, 10),
  mu = 2.0, u_fn = my_u_fn, u_max = 1.0,
  A = 0.2, alpha = 1.1, c_par = 0.02, p_par = 1.2,
  D_par = 0.05, gamma = 0.8, q_par = 1.5, beta = 2.0, m_min = 2.5
)

fun_rosen <- function(x, a = 1, b = 200) {
    return((a - x[1])^(2) + b * (x[2] - x[1]^(2))^(2))
}

grid_teste <- expand.grid(x = seq(-2, 2, length.out = 100),
                          y = seq(-2, 2, length.out = 100))

z <- matrix(apply(grid_teste, 1, fun_rosen),
            nrow = 100, ncol = 100)

image(x = seq(-2, 2, length.out = 100),
      y = seq(-2, 2, length.out = 100),
      z = z)


library(terra)

# 1. Obter a matriz do resultado do seu ajuste
matriz_bruta <- aj_peru$raster

# 2. Corrigir a orientação (transpor X/Y e inverter o eixo Y para Norte->Sul)
# Utiliza as variáveis y_grid e x_grid já definidas no seu ambiente
matriz_corrigida <- t(matriz_bruta)[length(y_grid):1, ]

# 3. Criar o SpatRaster definindo a extensão espacial correta
extensao <- terra::ext(min(x_grid), max(x_grid), min(y_grid), max(y_grid))
raster_bg <- terra::rast(matriz_corrigida, extent = extensao)

# OPCÃO A: Interpolação Visual (Mais rápida)
# O argumento smooth = TRUE força a interpolação bilinear na renderização do plot
terra::plot(raster_bg, smooth = TRUE,
            main = "Taxa de Fundo Contínua (Visual)")

# OPÇÃO B: Interpolação Física (Desagregação)
# Se precisar do mapa suave exportado como dados reais para outras análises,
# aumente a resolução do raster (ex: 4x mais pixéis) usando o método bilinear.
raster_alta_res <- terra::disagg(raster_bg, fact = 4, method = "bilinear")
terra::plot(raster_alta_res,
            main = "Taxa de Fundo Contínua (Dados Desagregados)")

# =============================================================================
# PATCH: background sem dependência da resolução da grade
#
# Problemas tratados
#  (1) d_min vs. espaçamento da grade: núcleos com h_j < Delta ficam mal
#      amostrados quando o KDE é avaliado apenas nos nós da grade.
#  (2) Hipótese  T * integral(mu_hat) = soma(rho_j): passa a valer EXATAMENTE
#      na janela (retângulo da grade), por construção.
#
# O que muda
#  - estimate_background_kernel_celula(): KDE por MASSA DE CÉLULA (diferenças de
#    pnorm nas bordas das células) + correção de borda c_j. Conserva a massa
#    para qualquer h_j, inclusive h_j << Delta. Vetorizado por produto de matrizes.
#  - mu_em_eventos(): mu_hat avaliado EXATAMENTE nos eventos (sem grade), com
#    opção leave-one-out (exclui a autocontribuição do evento j).
#  - calcular_probabilidades() e log_lik_etas(): aceitam o vetor mu_eventos
#    (e integral_bg), garantindo o MESMO mu nos passos E e M.
#  - ajustar_etas_proprio(): laço unificado usando as funções acima.
#
# Compatibilidade
#  - O raster devolvido continua na mesma grade e nas mesmas unidades
#    (eventos / unidade de tempo / unidade de área). Mapas, simulador e
#    resíduos continuam funcionando com ajuste$raster.
#  - Se mu_eventos = NULL, as funções usam func_nu() como antes.
#
# Convenção de grade: x_grid e y_grid são CENTROS de célula, com espaçamento
# constante; a janela é [x_1 - dx/2, x_N + dx/2] x [y_1 - dy/2, y_M + dy/2].
#
# Observação: este código foi escrito sem execução em R no ambiente de
# redação; a álgebra do estimador (conservação de massa, produto de matrizes)
# foi verificada numericamente em Python. Teste em um catálogo pequeno antes
# de rodar o catálogo completo.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. KDE adaptativo por massa de célula, com correção de borda
#    mu(x,y) = (1/T) sum_j (rho_j / c_j) K_{h_j}((x,y) - (x_j,y_j)),
#    c_j = massa do núcleo j dentro da janela.
#    Cada célula recebe a MASSA exata do núcleo (diferença de pnorm), dividida
#    pela área da célula => densidade média na célula.
#    Vale: sum(raster) * dx * dy == sum(rho_j) / T   (se corrigir_borda = TRUE)
# -----------------------------------------------------------------------------
estimate_background_kernel_celula <- function(x_grid, y_grid, x, y, rho_bg, bw_vec,
                                              T_total, corrigir_borda = TRUE,
                                              bloco = 2000) {
    Nx <- length(x_grid); Ny <- length(y_grid)
    dx <- diff(x_grid)[1]; dy <- diff(y_grid)[1]
    bx <- c(x_grid - dx / 2, x_grid[Nx] + dx / 2)   # bordas das células (Nx + 1)
    by <- c(y_grid - dy / 2, y_grid[Ny] + dy / 2)   # bordas das células (Ny + 1)
    n  <- length(x)
    # Massa de cada núcleo dentro da janela (núcleo gaussiano isotrópico: produto)
    cj <- (pnorm((bx[Nx + 1] - x) / bw_vec) - pnorm((bx[1] - x) / bw_vec)) *
          (pnorm((by[Ny + 1] - y) / bw_vec) - pnorm((by[1] - y) / bw_vec))
    cj <- pmax(cj, 1e-12)
    divisor <- if (corrigir_borda) cj else rep(1, n)
    w <- rho_bg / (T_total * divisor)               # peso de cada núcleo
    Z <- matrix(0, nrow = Nx, ncol = Ny)
    for (ini in seq(1, n, by = bloco)) {            # blocos limitam o uso de memória
        idx <- ini:min(ini + bloco - 1, n)
        m   <- length(idx)
        Px  <- pnorm(outer(bx, x[idx], "-") / rep(bw_vec[idx], each = Nx + 1))
        Py  <- pnorm(outer(by, y[idx], "-") / rep(bw_vec[idx], each = Ny + 1))
        Ax  <- Px[-1, , drop = FALSE] - Px[-(Nx + 1), , drop = FALSE]   # Nx x m
        Ay  <- Py[-1, , drop = FALSE] - Py[-(Ny + 1), , drop = FALSE]   # Ny x m
        Z   <- Z + (Ax * rep(w[idx], each = Nx)) %*% t(Ay)              # soma de produtos externos
    }
    list(raster = Z / (dx * dy),   # densidade (massa da célula / área)
         w = w,                    # pesos rho_j / (T c_j), reutilizados em mu_em_eventos()
         cj = cj)
}
# -----------------------------------------------------------------------------
# 2. mu_hat avaliado exatamente nos eventos (sem grade)
#    mu_j = sum_{l != j} w_l K_{h_l}((x_j,y_j) - (x_l,y_l))      (loo = TRUE)
#    Custo O(n^2) por iteração (paralelizável em j), da ordem do passo E.
#    O leave-one-out remove a autocontribuição w_j / (2 pi h_j^2), que
#    realimentaria rho_j na iteração seguinte quando h_j é pequeno.
# -----------------------------------------------------------------------------
mu_em_eventos <- function(x, y, w, bw_vec, loo = TRUE, num_cores = 1) {
    n      <- length(x)
    pre    <- w / (2 * pi * bw_vec^2)
    inv2h2 <- 1 / (2 * bw_vec^2)
    vals <- mclapply(seq_len(n), function(j) {
        d2 <- (x - x[j])^2 + (y - y[j])^2
        k  <- pre * exp(-d2 * inv2h2)
        if (loo) k[j] <- 0
        sum(k)
    }, mc.cores = num_cores)
    pmax(unlist(vals), 1e-15)                       # evita mu = 0 (log e razões)
}

# -----------------------------------------------------------------------------
# 3. Passo E (agora aceita mu_eventos)
# -----------------------------------------------------------------------------
calcular_probabilidades <- function(df_eventos, theta, M0, raster_mu, x_grid, y_grid,
                                    num_cores = 1, mu_eventos = NULL) {
    n <- nrow(df_eventos)
    mu_vals <- if (is.null(mu_eventos)) {
        func_nu(raster_mu, df_eventos, x_grid, y_grid)$lyr.1
    } else {
        mu_eventos
    }
    mu_vals <- pmax(mu_vals, 1e-15)
    probabilidades <- mclapply(1:n, function(j) {
        if (j == 1) return(numeric(0))
        i_indices <- 1:(j - 1)
        numeradores <- kappa_func(df_eventos$mag[i_indices], M0, theta["A"], theta["alpha"]) *
            g_func(df_eventos$time[j] - df_eventos$time[i_indices], theta["c"], theta["p"]) *
            f_spatial(
                (df_eventos$x[j] - df_eventos$x[i_indices])^2 + (df_eventos$y[j] - df_eventos$y[i_indices])^2,
                df_eventos$mag[i_indices], M0, theta["D"], theta["q"], theta["gamma"]
            )
        denom <- mu_vals[j] + sum(numeradores)
        return(numeradores / denom)
    }, mc.cores = num_cores)
    prob_total <- vapply(probabilidades, sum, numeric(1))
    list(probabilidades = probabilidades, prob_total = prob_total)
}

# -----------------------------------------------------------------------------
# 4. Log-verossimilhança (agora aceita mu_eventos e integral_bg)
#    integral_bg: valor de T * integral(mu). Com o KDE por célula e correção de
#    borda, sum(rho_j) é exatamente esse valor; o argumento permite passar outro.
# -----------------------------------------------------------------------------
log_lik_etas <- function(params, df_eventos, matriz_p_ij, prob_total, raster_mu, x_grid, y_grid,
                         M0 = 4.5, T_max = 10, mu_eventos = NULL, integral_bg = NULL) {
    n <- nrow(df_eventos)
    log_disparo_total <- 0
    mu_base_valores <- if (is.null(mu_eventos)) {
        func_nu(raster_mu, df_eventos, x_grid, y_grid)[, 1]
    } else {
        mu_eventos
    }
    for (j in 1:n) {
        rho_j <- pmax(0, 1 - prob_total[j])
        mu_j <- mu_base_valores[j]
        log_contrib_j <- rho_j * log(pmax(mu_j, 1e-15))
        if (j > 1) {
            p_ij <- matriz_p_ij[[j]]
            if (length(p_ij) > 0 && sum(p_ij) > 0) {
                i_indices <- 1:(j - 1)
                g_val <- g_func(df_eventos$time[j] - df_eventos$time[i_indices], params["c"], params["p"])
                kappa_val <- kappa_func(df_eventos$mag[i_indices], M0, params["A"], params["alpha"])
                f_val <- f_spatial(
                    (df_eventos$x[j] - df_eventos$x[i_indices])^2 + (df_eventos$y[j] - df_eventos$y[i_indices])^2,
                    df_eventos$mag[i_indices], M0, params["D"], params["q"], params["gamma"]
                )
                intensidade_ij <- kappa_val * g_val * f_val
                log_contrib_j <- log_contrib_j + sum(p_ij * log(pmax(intensidade_ij, 1e-15)))
            }
        }
        log_disparo_total <- log_disparo_total + log_contrib_j
    }
    integral_background <- if (is.null(integral_bg)) sum(pmax(0, 1 - prob_total)) else integral_bg
    integral_clustering <- sum(kappa_func(df_eventos$mag, M0, params["A"], params["alpha"]))
    log_lik_completa <- log_disparo_total - (integral_background + integral_clustering)
    return(-log_lik_completa)
}

# -----------------------------------------------------------------------------
# 5. Função principal de ajuste (laço unificado)
#    Novos argumentos:
#      corrigir_borda : KDE com correção de borda (massa exata na janela)
#      mu_exato       : mu nos eventos avaliado sem grade (senão, func_nu bilinear)
#      loo            : exclui o próprio evento ao avaliar mu nos eventos
#    Mudança de padrão: r0 constante = n / (T * area), em vez de matrix(1, ...)
#    (evita rho_j ~ 1 para todos os eventos no primeiro passo E).
# -----------------------------------------------------------------------------
ajustar_etas_proprio <- function(
                                 catalogo,
                                 theta0,
                                 M0 = 4.5,
                                 x_grid, y_grid,
                                 d_j_pre = NULL,
                                 r0 = NULL,
                                 num_cores = 1,
                                 max_iter_em = 100,
                                 tol_rel_theta  = 1e-3, tol_abs_theta  = 1e-6,
                                 tol_rel_raster = 1e-3, tol_abs_raster = 1e-6,
                                 lower = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-6),
                                 upper = c(15, 10, 15, 15, 15, 15, 10),
                                 guardar_historico = FALSE,
                                 corrigir_borda = TRUE,
                                 mu_exato = TRUE,
                                 loo = TRUE,
                                 retornar_pij = FALSE   # TRUE devolve a lista completa de p_ij (n^2/2 valores: pesada)
                                 ) {
    df_eventos <- extrair_eventos_catalog(catalogo, M0)
    df_eventos$time <- df_eventos$time - min(df_eventos$time)
    T_total <- max(df_eventos$time)
    n  <- nrow(df_eventos)
    dx <- diff(x_grid)[1]; dy <- diff(y_grid)[1]
    area <- length(x_grid) * dx * length(y_grid) * dy
    if (is.null(d_j_pre)) d_j_pre <- d_j(df_eventos[, c("x", "y")])
    if (is.null(r0)) r0 <- matrix(n / (T_total * area), nrow = length(x_grid), ncol = length(y_grid))
    # Um passo completo: E, KDE de mu, mu nos eventos, M
    passo_em <- function(theta_atual, r_atual, mu_ev_atual) {
        passo_e <- calcular_probabilidades(df_eventos, theta_atual, M0, r_atual, x_grid, y_grid,
                                           num_cores, mu_eventos = mu_ev_atual)
        rho_bg <- pmax(0, 1 - passo_e$prob_total)
        bg <- estimate_background_kernel_celula(x_grid, y_grid, df_eventos$x, df_eventos$y,
                                                rho_bg, d_j_pre, T_total, corrigir_borda)
        mu_ev_novo <- if (mu_exato) {
                          mu_em_eventos(df_eventos$x, df_eventos$y, bg$w, d_j_pre, loo, num_cores)
                      } else {
                          func_nu(bg$raster, df_eventos, x_grid, y_grid)$lyr.1
                      }
        otim <- optim(theta_atual, log_lik_etas, gr = gradiente_etas_completo,
                      df_eventos = df_eventos, matriz_p_ij = passo_e$probabilidades,
                      prob_total = passo_e$prob_total, raster_mu = bg$raster,
                      x_grid = x_grid, y_grid = y_grid, M0 = M0,
                      mu_eventos = mu_ev_novo, integral_bg = sum(rho_bg),
                      method = "L-BFGS-B", lower = lower, upper = upper,
                      control = list(maxit = 1000), hessian = TRUE)
        list(otim = otim, raster = bg$raster, mu_ev = mu_ev_novo)
    }
    thetas_hist  <- list(theta0)
    rasters_hist <- list(r0)
    theta_atual  <- theta0
    r_atual      <- r0
    mu_ev_atual  <- func_nu(r0, df_eventos, x_grid, y_grid)$lyr.1
    conv_codes   <- integer(0)
    i <- 1
    repeat {
        res <- passo_em(theta_atual, r_atual, mu_ev_atual)
        conv_codes     <- c(conv_codes, res$otim$convergence)
        theta_atual    <- res$otim$par
        hessiana_final <- res$otim$hessian
        r_atual        <- res$raster
        mu_ev_atual    <- res$mu_ev
        i <- i + 1
        thetas_hist[[i]]  <- theta_atual
        rasters_hist[[i]] <- r_atual
        cat("====Iteração", i - 1, "finalizada====\n")
        if (convergiu_em(rasters_hist[[i]], rasters_hist[[i - 1]], tol_rel_raster, tol_abs_raster) &&
            convergiu_em(thetas_hist[[i]], thetas_hist[[i - 1]], tol_rel_theta, tol_abs_theta)) {
            break
        }
        if (i >= max_iter_em) {
            warning(sprintf("ajustar_etas_proprio: max_iter_em = %d atingido sem convergir.", max_iter_em))
            break
        }
    }
    # E-step final em (theta_hat, mu_tilde final): probabilidade de cada evento ser de
    # background (rho_bg) e de ser disparado (prob_disparo = soma_i p_ij), na ordem de df_eventos.
    # Atenção ao nome: no código antigo, a coluna `rho` do catálogo era a prob. de ser DISPARADO;
    # aqui rho_bg é a de ser BACKGROUND (o rho_j do texto).
    pe_final <- calcular_probabilidades(df_eventos, theta_atual, M0, r_atual, x_grid, y_grid,
                                        num_cores, mu_eventos = mu_ev_atual)
    rho_bg <- pmax(0, 1 - pe_final$prob_total)
    pai_mais_provavel <- vapply(seq_len(n), function(j) {
        p <- pe_final$probabilidades[[j]]
        if (length(p) == 0 || rho_bg[j] >= max(p)) 0L else as.integer(which.max(p))
    }, integer(1))                       # 0 = background; i = índice do pai mais provável
    eventos <- cbind(df_eventos, rho_bg = rho_bg, prob_disparo = pe_final$prob_total,
                     pai_mais_provavel = pai_mais_provavel)
    list(
        theta              = theta_atual,
        raster             = r_atual,
        hessiana           = hessiana_final,
        mu_eventos         = mu_ev_atual,     # mu usado no passo E final (novo)
        rho_bg             = rho_bg,          # P(evento j é background), na ordem de `eventos` (novo)
        eventos            = eventos,         # df_eventos + rho_bg, prob_disparo, pai_mais_provavel (novo)
        probabilidades     = if (retornar_pij) pe_final$probabilidades else NULL,  # lista de p_ij (novo)
        bw                 = d_j_pre,         # bandas h_j (novo)
        T_total            = T_total,         # novo
        config             = list(corrigir_borda = corrigir_borda, mu_exato = mu_exato, loo = loo),
        n_iter_em          = i,
        convergence_codes  = conv_codes,
        convergiu          = all(conv_codes == 0) && i < max_iter_em,
        thetas_hist        = if (guardar_historico) thetas_hist  else NULL,
        rasters_hist       = if (guardar_historico) rasters_hist else NULL
    )
}

# -----------------------------------------------------------------------------
# 6. Diagnóstico do background (para o texto e para a análise de sensibilidade)
# -----------------------------------------------------------------------------
diagnosticar_background <- function(ajuste, df_eventos, x_grid, y_grid, M0 = 4.5, num_cores = 1) {
    dx <- diff(x_grid)[1]; dy <- diff(y_grid)[1]
    h  <- ajuste$bw
    pe <- calcular_probabilidades(df_eventos, ajuste$theta, M0, ajuste$raster, x_grid, y_grid,
                                  num_cores, mu_eventos = ajuste$mu_eventos)
    rho_bg <- pmax(0, 1 - pe$prob_total)
    mu_bilinear <- func_nu(ajuste$raster, df_eventos, x_grid, y_grid)$lyr.1
    list(
        espacamento           = c(dx = dx, dy = dy),
        frac_h_menor_meio_dx  = mean(h < max(dx, dy) / 2),   # núcleos mal resolvidos pela grade
        frac_h_no_piso        = mean(h <= min(h) + 1e-12),   # eventos no piso d_min
        razao_massa           = sum(ajuste$raster) * dx * dy / (sum(rho_bg) / ajuste$T_total),  # deve ser 1
        razao_mu_exato_bilin  = quantile(ajuste$mu_eventos / pmax(mu_bilinear, 1e-15),
                                         probs = c(0.05, 0.25, 0.5, 0.75, 0.95)),
        frac_background       = mean(rho_bg)
    )
}


# =============================================================================
# ALTERAÇÕES PONTUAIS nas funções de perfil (para usar o mesmo mu do ajuste)
# =============================================================================
# perfilar_rapido():
#   passo_e <- calcular_probabilidades(df_eventos, ajuste$theta, M0, ajuste$raster,
#                                      x_grid, y_grid, num_cores,
#                                      mu_eventos = ajuste$mu_eventos)
#   ... -log_lik_etas(theta_v, df_eventos, passo_e$probabilidades, passo_e$prob_total,
#                     ajuste$raster, x_grid, y_grid, M0 = M0,
#                     mu_eventos = ajuste$mu_eventos,
#                     integral_bg = sum(pmax(0, 1 - passo_e$prob_total)))
#
# perfilar_completo():
#   passo_e <- calcular_probabilidades(df_eventos, theta_ponto, M0, r_fixo,
#                                      x_grid, y_grid, num_cores_e,
#                                      mu_eventos = ajuste$mu_eventos)
#   obj <- function(theta_livre) { ...
#       log_lik_etas(theta_completo, df_eventos, passo_e$probabilidades, passo_e$prob_total,
#                    r_fixo, x_grid, y_grid, M0 = M0,
#                    mu_eventos = ajuste$mu_eventos,
#                    integral_bg = sum(pmax(0, 1 - passo_e$prob_total)))
#   }
#   (grad_obj não muda: gradiente_etas_completo não usa mu)
#
# Bootstrap paramétrico / resíduos: onde houver chamada a calcular_probabilidades()
# ou log_lik_etas() com ajuste$raster, passe também mu_eventos = ajuste$mu_eventos
# (ou, para catálogos simulados, o mu_eventos devolvido por ajustar_etas_proprio()).
# =============================================================================

load("~/resultados_sensibilidade_param.RData")

# 1. Obter a matriz do resultado do seu ajuste
matriz_correta <- raster_verd

# 2. Corrigir a orientação (transpor X/Y e inverter o eixo Y para Norte->Sul)
# Utiliza as variáveis y_grid e x_grid já definidas no seu ambiente
matriz_corrigida <- t(matriz_correta)[length(y_grid):1, ]

# 3. Criar o SpatRaster definindo a extensão espacial correta
extensao <- terra::ext(min(x_grid), max(x_grid), min(y_grid), max(y_grid))
raster_bg <- terra::rast(matriz_corrigida, extent = extensao)

# OPCÃO A: Interpolação Visual (Mais rápida)
# O argumento smooth = TRUE força a interpolação bilinear na renderização do plot
par(mfrow = c(2, 2))
terra::plot(raster_bg, smooth = TRUE,
            main = "True background rate")
terra::plot(t(resultados_sens[[1]]$raster)[length(y_grid):1, ] %>%  terra::rast(extent = extensao),
            smooth = T, main = "Reference background")
terra::plot(t(resultados_sens[[2]]$raster)[length(y_grid):1, ] %>%  terra::rast(extent = extensao),
            smooth = T, main = "Background of one LHS iteration")
terra::plot(t(resultados_sens[[126]]$raster)[length(y_grid):1, ] %>%  terra::rast(extent = extensao),
            smooth = T,  main = "Background of one MU0 iteration")
