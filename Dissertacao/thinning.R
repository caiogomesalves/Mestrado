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
theta_peru <- c("A" = 0.5, "alpha" = 1.0, "c" = 0.2,
                "p" = 1.15, "D" = 0.002, "beta" = 2)

# Simulação dos dados:
set.seed(5050)

sim_dados_peru <- simulate_hawkes_etas_zhuang(
  T_max = 10,
  X_lim = c(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long)),
  Y_lim = c(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat)),
  mu = 0.01, u_fn = rast_peru, u_max = max(im_peru),
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

# Extrator de valores de raster:
func_nu <- function(raster, dados) {
    terra::extract(raster, dados)
}

# Intensidade Total:
lambda_total <- function(t, long, lat, mag, M0, catalog, params, raster) {
    # Valor do background:
    mu_val <- func_nu(raster, catalog[, c("x", "y")])
    # Extração dos ancestrais:
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
    return(mu_val$lyr.1[nrow(clustering)] + clustering$total[nrow(clustering)])
}

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

# Raster inicial, considerando mu(x, y) = 1:
r <- terra::rast(nrows = 128, ncols = 128,
                 xmin = min(peru.quakes$longlat.coord$long),
                 xmax = max(peru.quakes$longlat.coord$long),
                 ymin = min(peru.quakes$longlat.coord$lat),
                 ymax = max(peru.quakes$longlat.coord$lat))

values(r) <- 1

rel_clust <- r

# Data.frame com os eventos para funções:
df_temp <- data.frame(x = catalogo_simulado_peru$longlat.coord$long,
                      y = catalogo_simulado_peru$longlat.coord$lat,
                      mag = catalogo_simulado_peru$revents[, "mm"] + 4,
                      time = catalogo_simulado_peru$revents[, "tt"])

# Distâncias para kernel adaptativo:
d_j_simulado <- d_j(catalogo_simulado_peru$longlat.coord[, c("long", "lat")])

# Cria listas com rasteres:
rasters <- list()

rast_rel_clust <- list()

rasters[[1]] <- r

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
        4, df_temp, theta_peru, r
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

# Estima o background rate:
values(r) <- estimate_background_kernel(seq(min(catalogo_simulado_peru$longlat.coord$long), max(catalogo_simulado_peru$longlat.coord$long), length.out = 128),
                                        seq(min(catalogo_simulado_peru$longlat.coord$lat), max(catalogo_simulado_peru$longlat.coord$lat), length.out = 128),
                                        data.frame(x = sim_dados_peru$x,
                                                   y = sim_dados_peru$y,
                                                   rho = pmax(0, 1 - prob_total)),
                                        d_j_simulado,
                                        T_total = max(catalogo_simulado_peru$rtperiod))

values(rel_clust) <- estimate_relative_clustering(seq(min(catalogo_simulado_peru$longlat.coord$long), max(catalogo_simulado_peru$longlat.coord$long), length.out = 128),
                                                  seq(min(catalogo_simulado_peru$longlat.coord$lat), max(catalogo_simulado_peru$longlat.coord$lat), length.out = 128),
                                                  data.frame(x = sim_dados_peru$x,
                                                             y = sim_dados_peru$y,
                                                             rho = prob_total),
                                                  d_j_simulado)$C

# Cria raster com valores estimados, e os inverte (não sei ainda o porquê):
r <- flip(r, direction = "vertical")
rel_clust <- flip(rel_clust, direction = "vertical")

rast_rel_clust[[1]] <- rel_clust

par(mfrow = c(1, 2))
plot(r)
plot(rel_clust)

#----Itera até a convergência----

rasters[[2]] <- r

i <- 2

# Considerando tolerância de 10e-6:
while (any(abs(values(rasters[[i]]) - values(rasters[[i - 1]])) > 10e-6)) {
    # Calcula pares de probabilidades:
    probabilidades <- mclapply(1:n_events, function(j) {
        if (j == 1) return(numeric(0))
        i_indices <- 1:(j - 1)
        denom <- lambda_total(
            df_temp$time[j], df_temp$x[j], df_temp$y[j], df_temp$mag[j],
            4, df_temp, theta_peru, r
        )
        numeradores <- kappa_func(df_temp$mag[i_indices], 4, theta_peru["A"], theta_peru["alpha"]) *
            g_func(df_temp$time[j] - df_temp$time[i_indices], theta_peru["c"], theta_peru["p"]) *
            f_spatial(
            (df_temp$x[j] - df_temp$x[i_indices])^2 + (df_temp$y[j] - df_temp$y[i_indices])^2,
            df_temp$mag[i_indices], 4, theta_peru["D"], theta_peru["alpha"]
            )
        return(numeradores / denom)
    }, mc.cores = num_cores)
    # Calcula probabilidades totais:
    prob_total <- probabilidades %>%
        lapply(sum) %>%
        do.call(c, .)
    # Estima background rate:
    values(r) <- estimate_background_kernel(seq(min(catalogo_simulado_peru$longlat.coord$long), max(catalogo_simulado_peru$longlat.coord$long), length.out = 128),
                                            seq(min(catalogo_simulado_peru$longlat.coord$lat), max(catalogo_simulado_peru$longlat.coord$lat), length.out = 128),
                                            data.frame(x = sim_dados_peru$x,
                                                       y = sim_dados_peru$y,
                                                       rho = prob_total),
                                            d_j_simulado,
                                            T_total = max(catalogo_simulado_peru$rtperiod))
    r <- flip(r, direction = "vertical")
    rasters[[i + 1]] <- r
    # Estima clustering relativo:
    values(rel_clust) <- estimate_relative_clustering(seq(min(catalogo_simulado_peru$longlat.coord$long), max(catalogo_simulado_peru$longlat.coord$long), length.out = 128),
                                                  seq(min(catalogo_simulado_peru$longlat.coord$lat), max(catalogo_simulado_peru$longlat.coord$lat), length.out = 128),
                                                  data.frame(x = sim_dados_peru$x,
                                                             y = sim_dados_peru$y,
                                                             rho = prob_total),
                                                  d_j_simulado)$C
    rel_clust <- flip(rel_clust, direction = "vertical")
    rast_rel_clust[[i]] <- rel_clust
    # Mostra progresso do loop:
    print(i + 1)
    # Atualiza índice:
    i <- i + 1
}

# Gráficos

par(mfrow = c(ceiling(sqrt(i)), ceiling(sqrt(i))))

# Background rate:
for (j in 1:i) {
    plot(rasters[[j]])
}

par(mfrow = c(ceiling(sqrt(i)), ceiling(sqrt(i - 1))))

# Clustering relativo:
for (j in 1:(i - 1)) {
    plot(rast_rel_clust[[j]])
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
