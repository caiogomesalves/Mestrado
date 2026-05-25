#----Pacotes necessários----

library(tidyverse)
library(ETAS)
library(spatstat)
library(readxl)
library(parallel)
library(terra)

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

#----Thinning----

# Parâmetros para estimação:
theta_peru <- c("A" = 0.5, "alpha" = 1.0, "c" = 0.2,
                "p" = 1.15, "D" = 0.002, "beta" = 2)

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
num_cores <- 8

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

print(paste("Tempo total de execução:", t2 - t1))

#----Salvando os Resultados----

save(
    rasters,
    rast_rel_clust,
    probabilidades,
    prob_total,
    t1, t2,
    file = "resultados_thinning_final.RData",
    compress = TRUE
)

if (exists("r") && !is.null(r)) {
    terra::writeRaster(r, "background_rate_final.tif", overwrite = TRUE)
}
if (exists("rel_clust") && !is.null(rel_clust)) {
    terra::writeRaster(rel_clust, "relative_clustering_final.tif", overwrite = TRUE)
}

print("Execução finalizada com sucesso. Arquivos salvos.")
