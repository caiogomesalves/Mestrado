#----Função K de Ripley----

func_k <- function(s, t, catalogo) {
    n <- nrow(catalogo$revents)
    range_x <- catalogo$region.win$xrange
    range_y <- catalogo$region.win$yrange
    A <- (range_x[2] - range_x[1]) * (range_y[2] - range_y[1])
    # D <- 1
    D <- catalogo$rtperiod[2]
    distancias <- as.matrix(dist(catalogo$longlat.coord[, c("lat", "long")]))
    tempos <- catalogo$revents[, "tt"]
    somas <- numeric(n)
    for (i in 1:n) {
        indices <- (1:n)[-i]
        tempos_diff <- abs(tempos[indices] - tempos[i])
        tempos_validos <- (tempos_diff < t)
        dist_validos <- (distancias[i, indices] < s)
        somas[i] <- sum((tempos_validos & dist_validos))
    }
    return((A * D/(n^2)) * sum(somas))
}

purrr::map(seq(0, 2, by = 0.01), function(x){func_k(x, 100, simu_peru)}) %>%
    do.call(c, .) %>%
    plot(seq(0, 2, by = 0.01), .)

lapply(seq(0.02, 0.5, by = 0.02), function(x){func_k(x, 10, simu_peru)})

func_k(0.145, 1, simu_peru)

Kest(as.ppp(simu_peru$revents[,c("xx","yy")],
            owin(simu_peru$region.win))) %>% plot()

#----Thinning----

d_j <- function(X, np = 20, dmin = 0.02) {
    distancias <- nndist(X, k = np)
    return(pmax(distancias, dmin))
}

d_j_peru <- d_j(simu_peru$longlat.coord[, c("long", "lat")])

# 1. Função de Resposta Espacial (f):
f_spatial <- function(dist, mag, M0, d, alpha) {
    sigma2 <- d * exp(alpha * (mag - M0))
    return((1 / (2 * pi * sigma2)) * exp(-(dist) / (2 * sigma2)))
}

func_nu <- function(raster, dados) {
    terra::extract(raster, dados)
}

# Intensidade Total (lambda_total):
lambda_total <- function(t, long, lat, mag, M0, catalog, params, raster) {
    mu_val <- func_nu(raster, catalog[, c("x", "y")])
    ancestors <- catalog %>% filter(time < t)
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

debugonce(lambda_total)

lambda_total(simu_peru$revents[, "tt"],
             simu_peru$longlat.coord[, "long"],
             simu_peru$longlat.coord[, "lat"],
             simu_peru$revents[, "mm"],
             4,
             catalog = data.frame(x = simu_peru$longlat.coord$long,
                                  y = simu_peru$longlat.coord$lat,
                                  time = simu_peru$revents[, "tt"]),
             theta_peru,
             rast_peru
             )$lyr.1

kappa_func(df_temp$mag[c(3537, 3538)], 4, 1, 1) *
    g_func(df_temp$time[3539] - df_temp$time[c(3537, 3538)], 0.2, 1.15) *
    f_spatial(
    (df_temp$x[3539] - df_temp$x[c(3537, 3538)])^2 + (df_temp$y[3539] - df_temp$y[c(3537, 3538)])^2,
    df_temp$mag[c(3537, 3538)], 4, 0.002, 1
    )

lambda_total(df_temp$time[3539],
             df_temp$x[3539],
             df_temp$y[3539],
             df_temp$mag[3539],
             4,
             df_temp,
             theta_peru2,
             r
             )

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

estimate_background_kernel(seq(min(simu_peru$longlat.coord$long), max(simu_peru$longlat.coord$long), length.out = 100),
                           seq(min(simu_peru$longlat.coord$lat), max(simu_peru$longlat.coord$lat), length.out = 100),
                           data.frame(x = simu_peru$longlat.coord$long, y = simu_peru$longlat.coord$lat, rho = rep(0.5, nrow(simu_peru$revents))),
                           d_j(catalogo_dados_homog$longlat.coord[, c("long", "lat")]),
                           T_total = 2000) %>%
    as.im() %>%
    plot()

probabilidades <- vector("list", nrow(catalogo_dados_homog_1$revents))

probabilidades[[1]] <- 0

# df_temp <- data.frame(x = cat_jap$longlat.coord$long,
#                       y = cat_jap$longlat.coord$lat,
#                       mag = cat_jap$revents[, "mm"] + 4,
#                       time = cat_jap$revents[, "tt"])

df_temp <- data.frame(x = catalogo_dados_homog_1$longlat.coord$long,
                      y = catalogo_dados_homog_1$longlat.coord$lat,
                      mag = catalogo_dados_homog_1$revents[, "mm"] + 4,
                      time = catalogo_dados_homog_1$revents[, "tt"])

theta_peru2 <- c("A" = 0.5, "alpha" = 1, "c" = 0.2,
                 "p" = 1.15, "D" = 0.002, "beta" = 2)

for (i in 1:nrow(catalogo_dados_homog_1$revents) - 1) {
    for (j in (i + 1):nrow(catalogo_dados_homog_1$revents)) {
        probabilidades[[j]][i] <- (kappa_func(df_temp$mag[i], 4, 2.5, 2) *
                                   g_func(df_temp$time[j] - df_temp$time[i], 1, 1.5) *
                                   f_spatial((df_temp$x[j] - df_temp$x[i])^2 + (df_temp$y[j] - df_temp$y[i])^2, df_temp$mag[i], 4,
                                             0.2, 0.1))/lambda_total(df_temp$time[j],
                                                                     df_temp$x[j],
                                                                     df_temp$y[j],
                                                                     df_temp$mag[j],
                                                                     4, df_temp, theta_peru, rast_peru)
    }
    print(i)
}

prob_total <- lapply(probabilidades, function(x) {sum(unlist(x))}) %>%
    do.call(c, .)

nu <- func_nu(rast_peru, data.frame(x = simu_peru$longlat.coord$long,
                                    y = simu_peru$longlat.coord$lat))$lyr.1

ggplot(sim_dados_peru, aes(x = x, y = y)) +
    geom_point(alpha = 0.1) +
    geom_point(data = sim_dados_peru %>% filter(parent == 374),
               aes(colour = gen))

#----Paralelização----

library(parallel)

# 1. Definir o número de núcleos (deixe 1 ou 2 livres para o sistema respirar)
num_cores <- max(1, detectCores() - 2)

# Número total de eventos:
n_events <- nrow(catalogo_dados_homog_1$revents)

# Paralelização do cálculo de probabilidades:
probabilidades_paralelas <- mclapply(1:n_events, function(j) {
    # Primeiro evento tem probabilidade 0 de ser descendente:
    if (j == 1) return(numeric(0))
    # Vetor de índices dos antecessores:
    i_indices <- 1:(j - 1)
    # Denominador comum para cada evento j:
    denom <- lambda_total(
        df_temp$time[j], df_temp$x[j], df_temp$y[j], df_temp$mag[j],
        4, df_temp, theta_peru2, r# rast_peru
    )
    # Computação do numerador para todos os antecessores:
    numeradores <- kappa_func(df_temp$mag[i_indices], 4, 0.5, 1) *
        g_func(df_temp$time[j] - df_temp$time[i_indices], 0.2, 1.15) *
        f_spatial(
        (df_temp$x[j] - df_temp$x[i_indices])^2 + (df_temp$y[j] - df_temp$y[i_indices])^2,
        df_temp$mag[i_indices], 4, 0.002, 1
        )
    # Retorna o vetor p_{i,j} para o evento j:
    return(numeradores / denom)
}, mc.cores = num_cores)

prob_total_par <- probabilidades_paralelas %>%
    lapply(sum) %>%
    do.call(c, .)

d_j_simu <- d_j(catalogo_dados_homog_1$longlat.coord[, c("long", "lat")])

values(r) <- estimate_background_kernel(seq(min(catalogo_dados_homog_1$longlat.coord$long), max(catalogo_dados_homog_1$longlat.coord$long), length.out = 128),
                                seq(min(catalogo_dados_homog_1$longlat.coord$lat), max(catalogo_dados_homog_1$longlat.coord$lat), length.out = 128),
                                data.frame(x = sim_dados_homog_1$x,
                                           y = sim_dados_homog_1$y,
                                           rho = prob_total_par),
                                d_j_simu,
                                T_total = 2000)

r <- flip(r, direction = "vertical")

plot(r)

#----Tamanho amostral----

calc_n <- function(alpha, beta, p1, p2) {
    pbar <- (p1 + p2)/2
    n = ((qnorm(1 - alpha/2) * sqrt(2 * pbar * (1 - pbar)) +
          qnorm(1 - beta) * sqrt(p1 * (1 - p1) + p2 * (1 - p2)))^2)/(p1 - p2)^2
    return(n)
}

calc_n(0.05, 0.1, 0.144, seq(0.15, 0.5, 0.05))

calc_n(0.05, 0.1, 0.144, 0.5)

par(mfrow = c(1, 1))

plot(seq(0.2, 0.6, 0.01),
     calc_n(0.01, 0.1, 0.144,
            seq(0.2, 0.6, 0.01)), ylim = c(0, 1000))

plot(seq(0.2, 0.6, 0.01),
     calc_n(0.05, 0.1, 0.144,
            seq(0.2, 0.6, 0.01)), ylim = c(0, 1000))

grid_n <- expand.grid(
    beta = seq(0.05, 0.25, 0.05),
    p_alt = seq(0.2, 0.95, 0.05)
)

grid_n %>%
    rowwise() %>%
    mutate(n = calc_n(0.05, beta, 0.144, p_alt) %>% round(0)) %>%
    ungroup() %>%
    pivot_wider(names_from = beta, values_from = n)

grid_n %>%
    rowwise() %>%
    mutate(n = calc_n(0.05, beta, 0.144, p_alt) %>% round(0)) %>%
    ungroup() %>%
    pivot_wider(names_from = beta, values_from = n)

grid_n %>%
    rowwise() %>%
    mutate(
        n = pwr.2p.test(h = ES.h(p_alt, 0.144), NULL, 0.05, 1 - beta, "greater")$n %>% round(0)
    ) %>%
    ungroup() %>%
    pivot_wider(names_from = beta, values_from = n)

for (i in seq(0.2, 0.95, 0.05)) {
    print(pwr.2p.test(h = ES.h(i, 0.144), n = NULL, sig.level = 0.05, power = 0.9, alternative = "greater")$n)
}

data.frame(
    x = seq(0, 1, 0.01)
) %>%
    rowwise() %>%
    mutate(
        y = ES.h(x, 0)
    ) %>%
    ungroup() %>%
    ggplot(aes(x = x, y = y)) +
    geom_line() +
    geom_segment(data = data.frame(x = 0.144, y = 0,
                                   xend = 0.144, yend = ES.h(0.144, 0),
                                   g = "Controle"),
                 aes(x = x, y = y, xend = xend, yend = yend,
                     colour = g)) +
    geom_segment(data = data.frame(x = 0.144, y = 0,
                                   xend = 0.144, yend = ES.h(0.144, 0),
                                   g = "Controle"),
                 aes(x = 0, y = yend, xend = xend, yend = yend,
                     colour = g)) +
    geom_segment(data = data.frame(x = 0.5, y = 0,
                                   xend = 0.5, yend = ES.h(0.5, 0),
                                   g = "Hemofílicos"),
                 aes(x = x, y = y, xend = xend, yend = yend,
                     colour = g)) +
    geom_segment(data = data.frame(x = 0.5, y = 0,
                                   xend = 0.5, yend = ES.h(0.5, 0),
                                   g = "Hemofílicos"),
                 aes(x = 0, y = yend, xend = xend, yend = yend,
                     colour = g)) +
    scale_x_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent,
                       limits = c(0, 1)) +
    labs(x = "Prevalências",
         y = "Prevalências transformadas",
         colour = "Grupo") +
    theme_bw() +
    theme(legend.position = "right")

pwr.2p.test(h = ES.h(0.5, 0.144), n = NULL, sig.level = 0.05, power = 0.9, alternative = "greater")

pwr.2p.test(h = NULL, n = 30, sig.level = 0.05, power = 0.9, alternative = "greater")

pwr.2p.test(h = ES.h(0.5, 0.144), n = NULL, sig.level = 0.05, power = 0.9, alternative = "greater") %>% plot()
