# Número total de eventos:
n_events <- nrow(peru.quakes$revents)

# Raster inicial, considerando mu(x, y) = 1:
r <- terra::rast(nrows = 128, ncols = 128,
                 xmin = min(peru.quakes$longlat.coord$long),
                 xmax = max(peru.quakes$longlat.coord$long),
                 ymin = min(peru.quakes$longlat.coord$lat),
                 ymax = max(peru.quakes$longlat.coord$lat))

values(r) <- 1

rel_clust <- r

# Data.frame com os eventos para funções:
df_temp <- data.frame(x = peru.quakes$longlat.coord$long,
                      y = peru.quakes$longlat.coord$lat,
                      mag = peru.quakes$revents[, "mm"] + 4,
                      time = peru.quakes$revents[, "tt"])

# Distâncias para kernel adaptativo:
d_j_peru <- d_j(peru.quakes$longlat.coord[, c("long", "lat")])

# Cria listas com rasteres:
rasters <- list()

rast_rel_clust <- list()

rasters[[1]] <- r

#----Estimação dos parents/offspring e de background rate----

# Definir o número de núcleos:
num_cores <- max(1, detectCores() - 4)

t1 <- Sys.time()

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
values(r) <- estimate_background_kernel(seq(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long), length.out = 128),
                                        seq(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat), length.out = 128),
                                        data.frame(x = peru$long,
                                                   y = peru$lat,
                                                   rho = pmax(0, 1 - prob_total)),
                                        d_j_peru,
                                        T_total = max(peru.quakes$rtperiod))

values(rel_clust) <- estimate_relative_clustering(seq(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long), length.out = 128),
                                                  seq(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat), length.out = 128),
                                                  data.frame(x = peru$long,
                                                             y = peru$lat,
                                                             rho = prob_total),
                                                  d_j_peru)$C

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
    values(r) <- estimate_background_kernel(seq(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long), length.out = 128),
                                            seq(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat), length.out = 128),
                                            data.frame(x = peru$long,
                                                       y = peru$lat,
                                                       rho = prob_total),
                                            d_j_peru,
                                            T_total = max(peru.quakes$rtperiod))
    r <- flip(r, direction = "vertical")
    rasters[[i + 1]] <- r
    # Estima clustering relativo:
    values(rel_clust) <- estimate_relative_clustering(seq(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long), length.out = 128),
                                                  seq(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat), length.out = 128),
                                                  data.frame(x = peru$long,
                                                             y = peru$lat,
                                                             rho = prob_total),
                                                  d_j_peru)$C
    rel_clust <- flip(rel_clust, direction = "vertical")
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
    plot(rasters[[j]])
}

par(mfrow = c(ceiling(sqrt(i)), ceiling(sqrt(i - 1))))

# Clustering relativo:
for (j in 1:(i - 1)) {
    plot(rast_rel_clust[[j]])
}

par(mfrow = c(1, 1))
