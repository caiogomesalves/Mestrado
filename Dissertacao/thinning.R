#----Pacotes necessários----

library(tidyverse)
library(ETAS)
library(spatstat)
library(readxl)
library(parallel)
library(terra)

#----Extração dos dados----

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

#----Simulação----

# Simulador baseado em Zhuang (2002):
simulate_hawkes_etas_zhuang <- function(
                                        # Parâmetros gerais:
                                        T_max = 100,           # Tempo máximo
                                        X_lim = c(-10, 10),    # Limites em X
                                        Y_lim = c(-10, 10),    # Limites em Y
                                        mu = 0.5,              # Taxa de background base
                                        u_fn = NULL,           # Função u(x,y) para o background não-homogêneo
                                        u_max = 1.0,           # Valor máximo de u_fn para aceitação/rejeição
                                        # Parâmetros ETAS:
                                        A = 0.1,               # Parâmetro de produtividade de escala
                                        alpha = 1.05,          # Eficiência de magnitude (Tempo e Espaço)
                                        c_par = 0.05,          # Deslocamento temporal (c)
                                        p_par = 1.2,           # Decaimento temporal de Omori (p > 1)
                                        d_par = 0.02,          # Escala de dispersão espacial (d)
                                        beta = 1.5,            # Parâmetro do Gutenberg-Richter
                                        m_min = 2.0            # Magnitude mínima (M0)
                                        ){
    # Lista para armazenar todos os eventos:
    all_events <- data.frame(t = numeric(), x = numeric(), y = numeric(), m = numeric(), gen = numeric())
    area <- (X_lim[2] - X_lim[1]) * (Y_lim[2] - Y_lim[1])
    # ---- 1. Geração dos Eventos Pais (Background / Imigrantes) ----
    # Caso seja dada uma função de background:
    if (!is.null(u_fn)) {
        lambda_max <- mu * u_max * area * T_max
        n_candidates <- rpois(1, lambda = lambda_max)
        if (n_candidates > 0) {
            c_x <- runif(n_candidates, X_lim[1], X_lim[2])
            c_y <- runif(n_candidates, Y_lim[1], Y_lim[2])
            mat_pontos <- matrix(c(c_x, c_y), ncol = 2)
            prob_accept <- terra::extract(u_fn, mat_pontos) / u_max
            prob_accept[is.na(prob_accept)] <- 0
            accepted <- runif(n_candidates) < prob_accept
            n_parents <- sum(accepted)
            if (n_parents > 0) {
                parents <- data.frame(
                    t = runif(n_parents, 0, T_max),
                    x = c_x[accepted],
                    y = c_y[accepted],
                    m = rexp(n_parents, rate = beta) + m_min,
                    gen = 0,
                    parent = 1:n_parents
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
        n_parents <- rpois(1, lambda = mu * area * T_max)
        if (n_parents > 0) {
            parents <- data.frame(
                t = runif(n_parents, 0, T_max),
                x = runif(n_parents, X_lim[1], X_lim[2]),
                y = runif(n_parents, Y_lim[1], Y_lim[2]),
                m = rexp(n_parents, rate = beta) + m_min,
                gen = 0,
                parent = 1:n_parents
            )
            all_events <- rbind(all_events, parents)
            current_gen <- parents
        } else {
            return(NULL)
        }
    }
    # ---- 2. Geração de Descendentes (Branching) ----
    generation <- 0
    while(nrow(current_gen) > 0) {
        generation <- generation + 1
        next_gen <- data.frame()
        for(i in 1:nrow(current_gen)) {
            parent <- current_gen[i,]
            # Produtividade: kappa(M) = A * exp(alpha * (M - M0)):
            expected_offspring <- A * exp(alpha * (parent$m - m_min))
            n_offspring <- rpois(1, expected_offspring)
            if(n_offspring > 0) {
                # Geração de Tempos (Inversa da CDF de Omori-Utsu):
                U_t <- runif(n_offspring)
                dt <- c_par * ((1 - U_t)^(-1 / (p_par - 1)) - 1)
                t_new <- parent$t + dt
                # Filtrar imediatamente os que passam do tempo máximo histórico:
                valid_t <- t_new <= T_max
                n_val <- sum(valid_t)
                if(n_val > 0) {
                    t_new <- t_new[valid_t]
                    # Geração do Espaço (Distribuição Gaussiana):
                    sigma2 <- d_par * exp(alpha * (parent$m - m_min))
                    sd_spatial <- sqrt(sigma2)
                    # Sorteio independente em X e Y simulando uma normal bivariada isotrópica:
                    x_new <- rnorm(n_val, mean = parent$x, sd = sd_spatial)
                    y_new <- rnorm(n_val, mean = parent$y, sd = sd_spatial)
                    # Marcas / Magnitudes
                    m_new <- rexp(n_val, rate = beta) + m_min
                    offspring <- data.frame(
                        t = t_new,
                        x = x_new,
                        y = y_new,
                        m = m_new,
                        gen = generation,
                        parent = parent$parent
                    )
                    # Filtro espacial para manter dentro da janela de interesse:
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
        if(generation > 100){
            warning("Regime supercrítico detectado (mais de 100 gerações). Interrompendo simulação.")
            break
        }
    }
    if(nrow(all_events) == 0) return(NULL)
    return(all_events[order(all_events$t), ])
}

# Raster baseado na base original:
im_peru <- ppp(x = peru$long, y = peru$lat,
               window = owin(xrange = c(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long)),
                             yrange = c(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat))),
               marks = peru$mag) %>%
    density.ppp(sigma = bw.scott)

rast_peru <- terra::rast(im_peru)

# Parâmetros para simulação:
theta_peru <- c("nu" = 0.01, "A" = 0.5, "alpha" = 1.0, "c" = 0.2,
                "p" = 1.15, "D" = 0.002, "beta" = 2)

# Simulação dos dados:
set.seed(5050)

sim_dados_peru <- simulate_hawkes_etas_zhuang(
  T_max = 10,
  X_lim = c(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long)),
  Y_lim = c(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat)),
  mu = theta_peru["nu"], u_fn = rast_peru, u_max = max(im_peru),
  A = theta_peru["A"], alpha = theta_peru["alpha"], c_par = theta_peru["c"],
  p_par = theta_peru["p"], d_par = theta_peru["D"], beta = theta_peru["beta"], m_min = 4
)

# Criação do catálogo:
catalogo_simulado_peru <- catalog(
    sim_dados_peru %>%
    mutate(t = dyears(t),
           data = POSIXct(1) + t,
           long = x, lat = y,
           mag = m,
           date = date(data),
           time = format(data, format = "%H:%M:%S"))
)

plot(catalogo_simulado_peru)

#----Funções auxiliares----

# Calculador de distâncias adaptativas:
d_j <- function(X, np = 20, dmin = 0.02) {
    distancias <- nndist(X, k = np)
    return(pmax(distancias, dmin))
}

# Produtividade:
kappa_func <- function(m, m0, A, alpha) {
    return(A * exp(alpha * (m - m0)))
}

# Densidade de probabilidade do tempo de ocorrência (Lei de Omori modificada):
g_func <- function(t_diff, c, p) {
    val <- rep(0, length(t_diff))
    idx <- t_diff > 0
    val[idx] <- (p - 1) * (c^(p - 1)) * (t_diff[idx] + c)^(-p)
    return(val)
}

# Densidade de resposta espacial:
f_spatial <- function(dist, mag, M0, d, alpha) {
    sigma2 <- d * exp(alpha * (mag - M0))
    return((1 / (2 * pi * sigma2)) * exp(-(dist) / (2 * sigma2)))
}

# # Extrator de valores de raster:
# func_nu <- function(raster, dados) {
#     terra::extract(raster, dados)
# }

# Extrator de valores de raster:
func_nu <- function(matriz_background, dados, x_grid, y_grid) {
    # Encontra o índice do grid mais próximo para cada ponto x e y:
    idx_x <- sapply(dados$x, function(val){which.min(abs(x_grid - val))})
    idx_y <- sapply(dados$y, function(val){which.min(abs(y_grid - val))})
    # Extrai os valores correspondentes da matriz:
    valores <- mapply(function(r, c) matriz_background[r, c], idx_x, idx_y)
    return(data.frame(lyr.1 = valores))
}

# Intensidade Total:
lambda_total <- function(t, long, lat, mag, M0, catalog, params, matriz_background, x_grid, y_grid) {
    # Valor do background:
    mu_val <- params["nu"] * func_nu(matriz_background, catalog[, c("x", "y")], x_grid, y_grid)
    ancestors <- catalog %>%
        filter(time < t)
    if(nrow(ancestors) == 0){
        return(mu_val)
    }
    clustering <- ancestors %>%
        mutate(
            dt = t - time,
            dist2 = (long - x)^2 + (lat - y)^2,
            contribution = kappa_func(mag, M0, params["A"], params["alpha"]) *
                g_func(dt, params["c"], params["p"]) *
                f_spatial(dist2, mag, M0, params["D"], params["alpha"])
        ) %>%
        summarise(total = sum(contribution))
    # Retorna o valor pontual da matriz + a soma do clustering:
    return(mu_val$lyr.1[nrow(clustering)] + clustering$total[nrow(clustering)])
}

# log_lik_etas <- function(params, df_eventos, matriz_p_ij, prob_total, raster_mu, x_grid, y_grid, M0 = 4, T_max = 10) {
#     n <- nrow(df_eventos)
#     nu_par <- params["nu"]
#     log_disparo_total <- 0
#     mu_base_valores <- func_nu(raster_mu, df_eventos, x_grid, y_grid)[, 1]
#     for (j in 1:n) {
#         rho_j <- pmax(0, 1 - prob_total[j])
#         mu_j <- nu_par * mu_base_valores[j]
#         log_contrib_j <- rho_j * log(pmax(mu_j, 1e-15))
#         if (j > 1) {
#             p_ij <- matriz_p_ij[[j]]
#             if (length(p_ij) > 0 && sum(p_ij) > 0) {
#                 i_indices <- 1:(j - 1)
#                 g_val <- g_func(df_eventos$time[j] - df_eventos$time[i_indices], params["c"], params["p"])
#                 kappa_val <- kappa_func(df_eventos$mag[i_indices], M0, params["A"], params["alpha"])
#                 f_val <- f_spatial(
#                 (df_eventos$x[j] - df_eventos$x[i_indices])^2 + (df_eventos$y[j] - df_eventos$y[i_indices])^2,
#                 df_eventos$mag[i_indices], M0, params["D"], params["alpha"]
#                 )
#                 intensidade_ij <- kappa_val * g_val * f_val
#                 log_contrib_j <- log_contrib_j + sum(p_ij * log(pmax(intensidade_ij, 1e-15)))
#             }
#         }
#         log_disparo_total <- log_disparo_total + log_contrib_j
#     }
#     integral_background <- nu_par * T_max
#     integral_clustering <- sum(kappa_func(df_eventos$mag, M0, params["A"], params["alpha"]))
#     log_lik_completa <- log_disparo_total - (integral_background + integral_clustering)
#     return(-log_lik_completa)
# }

log_lik_etas <- function(params, df_eventos, matriz_p_ij, prob_total, raster_mu, x_grid, y_grid, M0 = 4, T_max = 10) {
    n <- nrow(df_eventos)
    nu_par <- params["nu"]
    # 1. Termo dos Pares (Disparos Reais e Background)
    log_disparo_total <- 0
    # Extrai os valores do raster de background base para todos os eventos de uma vez
    mu_base_valores <- func_nu(raster_mu, df_eventos, x_grid, y_grid)[, 1]
    for (j in 1:n) {
        # --- Parte A: Contribuição do Background para o evento j ---
        rho_j <- pmax(0, 1 - prob_total[j])
        mu_j <- nu_par * mu_base_valores[j]
        log_contrib_j <- rho_j * log(pmax(mu_j, 1e-15))
        # --- Parte B: Contribuição do Clustering (para j > 1) ---
        if (j > 1) {
            p_ij <- matriz_p_ij[[j]]
            if (length(p_ij) > 0 && sum(p_ij) > 0) {
                i_indices <- 1:(j - 1)
                g_val <- g_func(df_eventos$time[j] - df_eventos$time[i_indices], params["c"], params["p"])
                kappa_val <- kappa_func(df_eventos$mag[i_indices], M0, params["A"], params["alpha"])
                f_val <- f_spatial(
                (df_eventos$x[j] - df_eventos$x[i_indices])^2 + (df_eventos$y[j] - df_eventos$y[i_indices])^2,
                df_eventos$mag[i_indices], M0, params["D"], params["alpha"]
                )
                intensidade_ij <- kappa_val * g_val * f_val
                log_contrib_j <- log_contrib_j + sum(p_ij * log(pmax(intensidade_ij, 1e-15)))
            }
        }
        log_disparo_total <- log_disparo_total + log_contrib_j
    }
    # Integral do background: nu * integral(mu(x,y)) * T_max
    integral_background <- nu_par * T_max
    # Integral do clustering:
    integral_clustering <- sum(kappa_func(df_eventos$mag, M0, params["A"], params["alpha"]))
    # Log-verossimilhança negativa
    log_lik_completa <- log_disparo_total - (integral_background + integral_clustering)
    return(-log_lik_completa)
}

gradiente_etas_completo <- function(params, df_eventos, matriz_p_ij, prob_total, T_max = 10, M0 = 4, ...) {
    nu_par <- params["nu"]
    A      <- params["A"]
    alpha  <- params["alpha"]
    cc     <- params["c"]
    pp     <- params["p"]
    D      <- params["D"]
    n <- nrow(df_eventos)
    # --- 1. Derivada em relação ao parâmetro nu ---
    # rho_j = 1 - prob_total_j
    rho <- pmax(0, 1 - prob_total)
    grad_nu <- -(sum(rho) / nu_par) + T_max
    # --- 2. Acumuladores para os parâmetros de Clustering (Passo E) ---
    sum_p_ij <- 0
    sum_alpha_term <- 0
    sum_c_term <- 0
    sum_p_term <- 0
    sum_D_term <- 0
    for (j in 2:n) {
        p_ij <- matriz_p_ij[[j]]
        if (length(p_ij) == 0 || sum(p_ij) == 0) next
        i_indices <- 1:(j - 1)
        dt <- df_eventos$time[j] - df_eventos$time[i_indices]
        dist2 <- (df_eventos$x[j] - df_eventos$x[i_indices])^2 + (df_eventos$y[j] - df_eventos$y[i_indices])^2
        dm <- df_eventos$mag[i_indices] - M0
        # Componente A
        sum_p_ij <- sum_p_ij + sum(p_ij)
        # Componente Alpha
        sigma2 <- D * exp(alpha * dm)
        deriv_alpha <- dm - (dist2 / (2 * sigma2))
        sum_alpha_term <- sum_alpha_term + sum(p_ij * deriv_alpha)
        # Componente c
        deriv_c <- ((pp - 1) / cc) - (pp / (dt + cc))
        sum_c_term <- sum_c_term + sum(p_ij * deriv_c)
        # Componente p
        deriv_p <- (1 / (pp - 1)) + log(cc) - log(dt + cc)
        sum_p_term <- sum_p_term + sum(p_ij * deriv_p)
        # Componente D
        deriv_D <- -(1 / D) + (dist2 / (2 * D * sigma2))
        sum_D_term <- sum_D_term + sum(p_ij * deriv_D)
    }
    # --- 3. Termos das Integrais de Clustering ---
    dm_all <- df_eventos$mag - M0
    exp_alpha_all <- exp(alpha * dm_all)
    int_A <- sum(exp_alpha_all)
    int_alpha <- sum(A * dm_all * exp_alpha_all)
    # --- 4. Montagem do Vetor do Gradiente Negativo ---
    grad_A     <- -(sum_p_ij / A) + int_A
    grad_alpha <- -sum_alpha_term + int_alpha
    grad_c     <- -sum_c_term
    grad_p     <- -sum_p_term
    grad_D     <- -sum_D_term
    return(c("nu" = grad_nu, "A" = grad_A, "alpha" = grad_alpha,
             "c" = grad_c, "p" = grad_p, "D" = grad_D))
}

# log_lik_etas <- function(params, df_eventos, matriz_p_ij, M0 = 4, T_max = 10) {
#     n <- nrow(df_eventos)
#     # 1. Termo dos Pares (Disparos Reais)
#     log_disparo_total <- 0
#     for (j in 2:n) {
#         p_ij <- matriz_p_ij[[j]] # Vetor de probabilidades para o evento j
#         if (length(p_ij) == 0 || sum(p_ij) == 0){
#             next
#         }
#         i_indices <- 1:(j - 1)
#         # Componentes do modelo usando os novos parâmetros testados pelo otimizador
#         g_val <- g_func(df_eventos$time[j] - df_eventos$time[i_indices], params["c"], params["p"])
#         kappa_val <- kappa_func(df_eventos$mag[i_indices], M0, params["A"], params["alpha"])
#         f_val <- f_spatial(
#         (df_eventos$x[j] - df_eventos$x[i_indices])^2 + (df_eventos$y[j] - df_eventos$y[i_indices])^2,
#         df_eventos$mag[i_indices], M0, params["D"], params["alpha"]
#         )
#         # Intensidade de clustering gerada de i para j
#         mu_ij <- kappa_val * g_val * f_val
#         # Evitar log(0) adicionando uma pequena tolerância
#         log_mu_ij <- log(pmax(mu_ij, 1e-15))
#         # Soma ponderada pelo peso p_ij atribuído no passo E
#         log_disparo_total <- log_disparo_total + sum(p_ij * log_mu_ij)
#     }
#     # 2. Termo Integral (Número esperado de descendentes que nasceram dentro do período)
#     # Nota: Para um modelo estrito, a integral de g_func sobre o tempo deveria ser considerada,
#     # mas uma aproximação comum assume a produtividade total kappa_func se T_max for grande.
#     integral_clustering <- sum(kappa_func(df_eventos$mag, M0, params["A"], params["alpha"]))
#     # Retorna a log-verossimilhança negativa (para minimizadores em R)
#     log_lik_completa <- log_disparo_total - integral_clustering
#     return(-log_lik_completa)
# }

# Estimação do background rate:
estimate_background_kernel <- function(x_grid, y_grid, catalog, bw_vec, T_total) {
    # Inicializa a matriz de densidade:
    z <- matrix(0, nrow = length(x_grid), ncol = length(y_grid))
    # Peso para cada evento: (1 - rho_j) / T:
    weights <- (1 - catalog$rho) / T_total
    # Loop sobre cada evento:
    for(i in 1:nrow(catalog)) {
        xi <- catalog$x[i]
        yi <- catalog$y[i]
        di <- bw_vec[i]
        wi <- weights[i]
        dx2 <- (x_grid - xi)^2
        dy2 <- (y_grid - yi)^2
        contribution <- (wi / (2 * pi * di^2)) * exp(-outer(dx2, dy2, "+") / (2 * di^2))
        z <- z + contribution
    }
    return(z)
}

# Estimação do efeito de clustering relativo:
estimate_relative_clustering <- function(x_grid, y_grid, catalog, bw_vec) {
    # Inicializa as matrizes para o numerador e denominador:
    numerator_gamma <- matrix(0, nrow = length(x_grid), ncol = length(y_grid))
    denominator_m1  <- matrix(0, nrow = length(x_grid), ncol = length(y_grid))
    # Loop sobre cada evento para acumular as superfícies:
    for(i in 1:nrow(catalog)) {
        xi <- catalog$x[i]
        yi <- catalog$y[i]
        di <- bw_vec[i]
        rho_i <- catalog$rho[i]
        dx2 <- (x_grid - xi)^2
        dy2 <- (y_grid - yi)^2
        # Kernel Gaussiano base avaliado no grid para o evento i:
        kernel_base <- (1 / (2 * pi * di^2)) * exp(-outer(dx2, dy2, "+") / (2 * di^2))
        # Acumula o denominador:
        denominator_m1 <- denominator_m1 + kernel_base
        # Acumula o numerador:
        numerator_gamma <- numerator_gamma + (rho_i * kernel_base)
    }
    # Calcula a razão ponto a ponto para obter C(x,y)
    C_xy <- ifelse(denominator_m1 > 0, numerator_gamma / denominator_m1, 0)
    return(list(C = C_xy, m1 = denominator_m1, gamma = numerator_gamma))
}

#----Algoritmo 1 do artigo Zhuang (2002)----

# Número total de eventos:
n_events <- nrow(catalogo_simulado_peru$revents)

# Definir as sequências do grid explicitamente
x_grid <- seq(min(catalogo_simulado_peru$longlat.coord$long), max(catalogo_simulado_peru$longlat.coord$long), length.out = 128)
y_grid <- seq(min(catalogo_simulado_peru$longlat.coord$lat), max(catalogo_simulado_peru$longlat.coord$lat), length.out = 128)

r <- matrix(1, nrow = 128, ncol = 128)
rel_clust <- matrix(0, nrow = 128, ncol = 128)

r <- r/(sum(r) * (xres(rast(r)) * yres(rast(r))))

# Na lista de rasters, guardaremos as matrizes:
rasters <- list()
rast_rel_clust <- list()
rasters[[1]] <- r

thetas <- list()
thetas[[1]] <- theta_peru

# # Raster inicial, considerando mu(x, y) = 1:
# r <- terra::rast(nrows = 128, ncols = 128,
#                  xmin = min(peru.quakes$longlat.coord$long),
#                  xmax = max(peru.quakes$longlat.coord$long),
#                  ymin = min(peru.quakes$longlat.coord$lat),
#                  ymax = max(peru.quakes$longlat.coord$lat))
#
# values(r) <- 1
#
# rel_clust <- r

# Data.frame com os eventos para funções:
df_temp <- data.frame(x = catalogo_simulado_peru$longlat.coord$long,
                      y = catalogo_simulado_peru$longlat.coord$lat,
                      mag = catalogo_simulado_peru$revents[, "mm"] + 4,
                      time = catalogo_simulado_peru$revents[, "tt"])

# Distâncias para kernel adaptativo:
d_j_simulado <- d_j(catalogo_simulado_peru$longlat.coord[, c("long", "lat")])

# Cria listas com rasteres:
# rasters <- list()
#
# rast_rel_clust <- list()
#
# rasters[[1]] <- r

#----Estimação dos parents/offspring e de background rate----

# Definir o número de núcleos:
num_cores <- max(1, detectCores() - 2)

# Paralelização do cálculo de probabilidades:
probabilidades <- mclapply(1:n_events, function(j) {
    # Primeiro evento tem probabilidade 0 de ser descendente:
    if (j == 1) return(numeric(0))
    # Vetor de índices dos antecessores:
    i_indices <- 1:(j - 1)
    # Denominador comum para cada evento j:
    denom <- lambda_total(
        df_temp$time[j], df_temp$x[j], df_temp$y[j], df_temp$mag[j],
        4, df_temp, theta_peru, r, x_grid, y_grid
    )
    # Computação do numerador para todos os antecessores:
    numeradores <- kappa_func(df_temp$mag[i_indices], 4, theta_peru["A"], theta_peru["alpha"]) *
        g_func(df_temp$time[j] - df_temp$time[i_indices], theta_peru["c"], theta_peru["p"]) *
        f_spatial(
        (df_temp$x[j] - df_temp$x[i_indices])^2 + (df_temp$y[j] - df_temp$y[i_indices])^2,
        df_temp$mag[i_indices], 4, theta_peru["D"], theta_peru["alpha"]
        )
    # Retorna o vetor p_{i,j} para o evento j:
    return(numeradores / denom)
}, mc.cores = num_cores)

# Probabilidades totais:
prob_total <- probabilidades %>%
    lapply(sum) %>%
    do.call(c, .)

otimizacao <- optim(theta_peru[names(theta_peru) != "beta"],
                    log_lik_etas,
                    gr = gradiente_etas_completo,
                    df_eventos = df_temp,
                    matriz_p_ij = probabilidades,
                    prob_total = prob_total,
                    raster_mu = r, x_grid = x_grid, y_grid = y_grid,
                    T_max = max(df_temp$time),
                    method = "L-BFGS-B",
                    lower = c(0.001, 0.001, 0.01, 0.001, 1.001, 0.0001),
                    upper = c(1e6, 10, 10, 10, 10, 10),
                    control = list(maxit = 1000))

theta_peru <- otimizacao$par
thetas[[2]] <- theta_peru

# Estima o background rate:
r <- estimate_background_kernel(seq(min(catalogo_simulado_peru$longlat.coord$long), max(catalogo_simulado_peru$longlat.coord$long), length.out = 128),
                                        seq(min(catalogo_simulado_peru$longlat.coord$lat), max(catalogo_simulado_peru$longlat.coord$lat), length.out = 128),
                                        data.frame(x = sim_dados_peru$x,
                                                   y = sim_dados_peru$y,
                                                   rho = pmax(0, 1 - prob_total)),
                                        d_j_simulado,
                                        T_total = max(catalogo_simulado_peru$rtperiod))

rel_clust <- estimate_relative_clustering(seq(min(catalogo_simulado_peru$longlat.coord$long), max(catalogo_simulado_peru$longlat.coord$long), length.out = 128),
                                                  seq(min(catalogo_simulado_peru$longlat.coord$lat), max(catalogo_simulado_peru$longlat.coord$lat), length.out = 128),
                                                  data.frame(x = sim_dados_peru$x,
                                                             y = sim_dados_peru$y,
                                                             rho = prob_total),
                                                  d_j_simulado)$C

# Cria raster com valores estimados, e os inverte (não sei ainda o porquê):
# r <- flip(r, direction = "vertical")
# rel_clust <- flip(rel_clust, direction = "vertical")

rast_rel_clust[[1]] <- rel_clust

par(mfrow = c(1, 2))
plot(as.im(r))
plot(as.im(rel_clust))

#----Itera até a convergência----

rasters[[2]] <- r

i <- 2

t1 <- Sys.time()

# Considerando tolerância de 10e-6:
while (any(abs(rasters[[i]] - rasters[[i - 1]]) > 10e-5)) {
    # Calcula pares de probabilidades:
    # Paralelização do cálculo de probabilidades:
    probabilidades <- mclapply(1:n_events, function(j) {
        # Primeiro evento tem probabilidade 0 de ser descendente:
        if (j == 1) return(numeric(0))
        # Vetor de índices dos antecessores:
        i_indices <- 1:(j - 1)
        # Denominador comum para cada evento j:
        denom <- lambda_total(
            df_temp$time[j], df_temp$x[j], df_temp$y[j], df_temp$mag[j],
            4, df_temp, theta_peru, r, x_grid, y_grid
        )
        # Computação do numerador para todos os antecessores:
        numeradores <- kappa_func(df_temp$mag[i_indices], 4, theta_peru["A"], theta_peru["alpha"]) *
            g_func(df_temp$time[j] - df_temp$time[i_indices], theta_peru["c"], theta_peru["p"]) *
            f_spatial(
            (df_temp$x[j] - df_temp$x[i_indices])^2 + (df_temp$y[j] - df_temp$y[i_indices])^2,
            df_temp$mag[i_indices], 4, theta_peru["D"], theta_peru["alpha"]
            )
        # Retorna o vetor p_{i,j} para o evento j:
        return(numeradores / denom)
    }, mc.cores = num_cores)
    # probabilidades <- mclapply(1:n_events, function(j) {
    #     if (j == 1) return(numeric(0))
    #     i_indices <- 1:(j - 1)
    #     denom <- lambda_total(
    #         df_temp$time[j], df_temp$x[j], df_temp$y[j], df_temp$mag[j],
    #         4, df_temp, theta_peru, r, x_grid, y_grid
    #     )
    #     numeradores <- kappa_func(df_temp$mag[i_indices], 4, theta_peru["A"], theta_peru["alpha"]) *
    #         g_func(df_temp$time[j] - df_temp$time[i_indices], theta_peru["c"], theta_peru["p"]) *
    #         f_spatial(
    #         (df_temp$x[j] - df_temp$x[i_indices])^2 + (df_temp$y[j] - df_temp$y[i_indices])^2,
    #         df_temp$mag[i_indices], 4, theta_peru["D"], theta_peru["alpha"]
    #         )
    #     return(numeradores / denom)
    # }, mc.cores = num_cores)
    # Calcula probabilidades totais:
    prob_total <- probabilidades %>%
        lapply(sum) %>%
        do.call(c, .)
    # Estima background rate:
    otimizacao <- optim(theta_peru[names(theta_peru) != "beta"],
                        log_lik_etas,
                        gr = gradiente_etas_completo,
                        df_eventos = df_temp,
                        matriz_p_ij = probabilidades,
                        prob_total = prob_total,
                        raster_mu = r, x_grid = x_grid, y_grid = y_grid,
                        T_max = max(df_temp$time),
                        method = "L-BFGS-B",
                        lower = c(0.001, 0.001, 0.01, 0.001, 1.001, 0.0001),
                        upper = c(1e6, 2, 2, 5, 5, 2),
                        control = list(maxit = 1000))
    theta_peru <- otimizacao$par
    thetas[[i + 1]] <- theta_peru
    r <- estimate_background_kernel(seq(min(catalogo_simulado_peru$longlat.coord$long), max(catalogo_simulado_peru$longlat.coord$long), length.out = 128),
                                    seq(min(catalogo_simulado_peru$longlat.coord$lat), max(catalogo_simulado_peru$longlat.coord$lat), length.out = 128),
                                    data.frame(x = sim_dados_peru$x,
                                               y = sim_dados_peru$y,
                                               rho = prob_total),
                                    d_j_simulado,
                                    T_total = max(catalogo_simulado_peru$rtperiod))
    # r <- flip(r, direction = "vertical")
    rasters[[i + 1]] <- r
    # Estima clustering relativo:
    rel_clust <- estimate_relative_clustering(seq(min(catalogo_simulado_peru$longlat.coord$long), max(catalogo_simulado_peru$longlat.coord$long), length.out = 128),
                                              seq(min(catalogo_simulado_peru$longlat.coord$lat), max(catalogo_simulado_peru$longlat.coord$lat), length.out = 128),
                                              data.frame(x = sim_dados_peru$x,
                                                         y = sim_dados_peru$y,
                                                         rho = prob_total),
                                              d_j_simulado)$C
    # rel_clust <- flip(rel_clust, direction = "vertical")
    rast_rel_clust[[i]] <- rel_clust
    # Mostra progresso do loop:
    print(i + 1)
    # Atualiza índice:
    i <- i + 1
}

t2 <- Sys.time()

# Gráficos

par(mfrow = c(ceiling(sqrt(i)), ceiling(sqrt(i))))

# Background rate:
for (j in 1:i) {
    plot(rasters[[j]] %>% as.im())
}

par(mfrow = c(ceiling(sqrt(i)), ceiling(sqrt(i - 1))))

# Clustering relativo:
for (j in 1:(i - 1)) {
    plot(rast_rel_clust[[j]] %>% as.im())
}

par(mfrow = c(1, 1))

#----Thinning----

set.seed(5050)
# Gera valores de uma uniforme:
U <- runif(n_events)

# Avalia quais foram mantidos, como pais:
sim_dados_peru[U < (1 - prob_total), ] %>%
    mutate(Pais = case_when(gen == 0 ~ "Pai",
                            gen != 0 ~ "Filho")) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(colour = Pais)) +
    scale_colour_manual(values = c("red", "steelblue")) +
    scale_x_continuous(limits = range(sim_dados_peru$x)) +
    scale_y_continuous(limits = range(sim_dados_peru$y)) +
    theme_bw()

# Avalia quais foram mantidos, como pais:
sim_dados_peru[U < (1 - prob_total), ] %>%
    mutate(Pais = case_when(gen == 0 ~ "Pai",
                            gen != 0 ~ "Filho")) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(colour = Pais)) +
    scale_colour_manual(values = c("red", "steelblue")) +
    scale_x_continuous(limits = range(sim_dados_peru$x)) +
    scale_y_continuous(limits = range(sim_dados_peru$y)) +
    theme_bw()

# Avalia quais foram mantidos, como filhos:
sim_dados_peru[U > (1 - prob_total), ] %>%
    mutate(Pais = case_when(gen == 0 ~ "Pai",
                            gen != 0 ~ "Filho")) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(colour = Pais)) +
    scale_colour_manual(values = c("red", "steelblue")) +
    scale_x_continuous(limits = range(sim_dados_peru$x)) +
    scale_y_continuous(limits = range(sim_dados_peru$y)) +
    theme_bw()

# Quantidade total de pais e filhos na base original:
sim_dados_peru %>%
    summarise(Pais = sum(gen == 0),
              Filhos = sum(gen != 0))

# Quantidade relativa de pais e filhos após thinning:
(sim_dados_peru[U < (1 - prob_total), ] %>%
 summarise(Pais = sum(gen == 0),
           Filhos = sum(gen != 0)))/(sim_dados_peru %>%
                                     summarise(Pais = sum(gen == 0),
                                               Filhos = sum(gen != 0)))


#----Teste----

teste <- rast(r, crs = "+proj=longlat +datum=WGS84")
teste <- flip(teste)
ext(teste) <- c(min(peru.quakes$longlat.coord$long),
                max(peru.quakes$longlat.coord$long),
                min(peru.quakes$longlat.coord$lat),
                max(peru.quakes$longlat.coord$lat))

# Simulação dos dados:
set.seed(5050)

sim_dados_peru_2 <- simulate_hawkes_etas_zhuang(
  T_max = 10,
  X_lim = c(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long)),
  Y_lim = c(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat)),
  mu = 150, u_fn = teste, u_max = max(r),
  A = 1, alpha = 1, c_par = theta_peru["c"],
  p_par = theta_peru["p"], d_par = theta_peru["D"], beta = theta_peru["beta"], m_min = 4
)

# Criação do catálogo:
catalogo_simulado_peru_2 <- catalog(
    sim_dados_peru_2 %>%
    mutate(t = dyears(t),
           data = POSIXct(1) + t,
           long = x, lat = y,
           mag = m,
           date = date(data),
           time = format(data, format = "%H:%M:%S"))
)

plot(catalogo_simulado_peru_2)

#----Resultados Cluster----

length(rasters)

par(mfrow = c(4, 4))

for (i in 1:length(rasters)) {
    plot(as.im(rasters[[i]]))
}

for (i in 1:(length(rasters) - 1)) {
    plot(as.im(rast_rel_clust[[i]]))
}

par(mfrow = c(1, 1))

hist(prob_total)

thetas %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    select(-c("nu")) %>%
    matplot(type = 'l')

set.seed(1)

U <- runif(nrow(peru))

pais <- peru[(U < 1 - prob_total), ]

filhos <- peru[(U >= 1 - prob_total), ]

pais %>%
    ggplot(aes(x = lat, y = long)) +
    geom_point()

filhos %>%
    ggplot(aes(x = lat, y = long)) +
    geom_point()

pais %>%
    catalog() %>%
    plot()

filhos %>%
    catalog() %>%
    plot()

#----Resíduos espaciais----

rast_teste <- (as.im(rasters[[14]],
                     owin(xrange = c(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long)),
                          yrange = c(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat)))) %>% rast())

sum_r <- sum(values(rast_teste))

normalizador <- sum_r * (res(rast_teste)[1] * res(rast_teste)[2])

# 1. Transformar seu raster de background (estimado no EM) em uma imagem do spatstat (im)
# Multiplicamos pelo Tempo Total porque 'r' é uma taxa temporal e queremos o valor esperado acumulado
mu_esperado <- as.im(rasters[[14]]/normalizador,
                     owin(xrange = c(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long)),
                          yrange = c(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat)))) * peru.quakes$rtperiod[2]

# 2. Criar um padrão de pontos (ppp) do catálogo com "marcas" sendo as probabilidades
# Aqui, prob_total é o (1 - rho) do seu EM, a probabilidade do evento i ser background
pontos_obs <- ppp(
    x = peru$long,
    y = peru$lat,
    window = owin(xrange = c(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long)),
                  yrange = c(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat))),
    marks = prob_total
)

# 3. Calcular a densidade espacial observada do background
# Usamos as probabilidades como "pesos" (weights) na suavização
mu_observado <- density.ppp(pontos_obs,
                            weights = pontos_obs$marks,
                            sigma = bw.diggle(pontos_obs))

# 4. Calcular o Raster de Resíduos
residuos_espaciais <- eval.im((mu_observado - mu_esperado)/sqrt(mu_esperado + 1e-10))

# 5. Plotar
plot(residuos_espaciais, main = "Mapa de Resíduos Espaciais do Background")
contour(residuos_espaciais, add = TRUE, col = "black", alpha = 0.5)
points(pontos_obs, pch = 20, cex = 0.2, col = rgb(0,0,0,0.3))
