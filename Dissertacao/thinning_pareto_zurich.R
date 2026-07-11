#----Pacotes necessários----

library(dplyr)
library(ETAS)
library(spatstat)
library(readxl)
library(parallel)
library(httr)

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

# Densidade de resposta espacial (Pareto/Lei de Potência - Pacote ETAS):
f_spatial <- function(dist2, mag, M0, D_par, q_par, gamma_par) {
    # S_i é o fator de escala espacial dependente da magnitude
    S_i <- D_par * exp(gamma_par * (mag - M0))
    # f(x,y) = ((q-1) / (pi * S_i)) * (1 + dist2 / S_i)^(-q)
    term1 <- (q_par - 1) / (pi * S_i)
    term2 <- (1 + (dist2 / S_i))^(-q_par)
    return(term1 * term2)
}

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
    # mu_val <- func_nu(matriz_background, catalog[, c("x", "y")], x_grid, y_grid)
    mu_val <- func_nu(matriz_background, data.frame(x = long, y = lat), x_grid, y_grid)
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
                f_spatial(dist2, mag, M0, params["D"], params["q"], params["gamma"])
        ) %>%
        summarise(total = sum(contribution))
    return(mu_val$lyr.1[nrow(clustering)] + clustering$total[nrow(clustering)])
}

# Calcula as probabilidades p_ij (passo E do EM) para todos os eventos de uma vez:
calcular_probabilidades <- function(df_eventos, theta, M0, raster_mu, x_grid, y_grid, num_cores = 1) {
    n <- nrow(df_eventos)
    # Valores de fundo (mu) para todos os eventos, calculados de uma só vez:
    mu_vals <- func_nu(raster_mu, df_eventos, x_grid, y_grid)$lyr.1
    probabilidades <- mclapply(1:n, function(j) {
        # Primeiro evento tem probabilidade 0 de ser descendente:
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

# Log-verossimilhança do modelo ajustada:
log_lik_etas <- function(params, df_eventos, matriz_p_ij, prob_total, raster_mu, x_grid, y_grid, M0 = 4.5, T_max = 10) {
    n <- nrow(df_eventos)
    log_disparo_total <- 0
    mu_base_valores <- func_nu(raster_mu, df_eventos, x_grid, y_grid)[, 1]
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
                # Chamada atualizada para o Kernel de Pareto
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
    integral_background <- sum(pmax(0, 1 - prob_total))
    integral_clustering <- sum(kappa_func(df_eventos$mag, M0, params["A"], params["alpha"]))
    log_lik_completa <- log_disparo_total - (integral_background + integral_clustering)
    return(-log_lik_completa)
}

# Gradiente Analítico Exato (incluindo q e gamma):
gradiente_etas_completo <- function(params, df_eventos, matriz_p_ij, prob_total, T_max = 10, M0 = 4.5, ...) {
    A      <- params["A"]
    alpha  <- params["alpha"]
    cc     <- params["c"]
    pp     <- params["p"]
    D      <- params["D"]
    qq     <- params["q"]
    gamma  <- params["gamma"]
    n <- nrow(df_eventos)
    sum_p_ij <- 0
    sum_alpha_term <- 0
    sum_c_term <- 0
    sum_p_term <- 0
    sum_D_term <- 0
    sum_q_term <- 0
    sum_gamma_term <- 0
    for (j in 2:n) {
        p_ij <- matriz_p_ij[[j]]
        if (length(p_ij) == 0 || sum(p_ij) == 0) next
        i_indices <- 1:(j - 1)
        dt <- df_eventos$time[j] - df_eventos$time[i_indices]
        dist2 <- (df_eventos$x[j] - df_eventos$x[i_indices])^2 + (df_eventos$y[j] - df_eventos$y[i_indices])^2
        dm <- df_eventos$mag[i_indices] - M0
        sum_p_ij <- sum_p_ij + sum(p_ij)
        # Derivada de Produtividade (alpha)
        sum_alpha_term <- sum_alpha_term + sum(p_ij * dm)
        # Derivadas Temporais
        deriv_c <- ((pp - 1) / cc) - (pp / (dt + cc))
        sum_c_term <- sum_c_term + sum(p_ij * deriv_c)
        deriv_p <- (1 / (pp - 1)) + log(cc) - log(dt + cc)
        sum_p_term <- sum_p_term + sum(p_ij * deriv_p)
        # Derivadas Espaciais (Kernel de Pareto)
        S_i <- D * exp(gamma * dm)
        V_ij <- dist2 / S_i
        U_ij <- 1 + V_ij
        # Fator de escala comum das derivadas em D e gamma
        W_ij <- qq * (V_ij / U_ij) - 1
        deriv_D <- W_ij / D
        sum_D_term <- sum_D_term + sum(p_ij * deriv_D)
        deriv_gamma <- dm * W_ij
        sum_gamma_term <- sum_gamma_term + sum(p_ij * deriv_gamma)
        deriv_q <- (1 / (qq - 1)) - log(U_ij)
        sum_q_term <- sum_q_term + sum(p_ij * deriv_q)
    }
    dm_all <- df_eventos$mag - M0
    exp_alpha_all <- exp(alpha * dm_all)
    int_A <- sum(exp_alpha_all)
    int_alpha <- sum(A * dm_all * exp_alpha_all)
    grad_A     <- -(sum_p_ij / A) + int_A
    grad_alpha <- -sum_alpha_term + int_alpha
    grad_c     <- -sum_c_term
    grad_p     <- -sum_p_term
    grad_D     <- -sum_D_term
    grad_q     <- -sum_q_term
    grad_gamma <- -sum_gamma_term
    return(c("A" = grad_A, "alpha" = grad_alpha,
             "c" = grad_c, "p" = grad_p, "D" = grad_D,
             "q" = grad_q, "gamma" = grad_gamma))
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

# Critério de convergência do EM:
convergiu_em <- function(novo, antigo, tol_rel = 1e-3, tol_abs = 1e-6) {
    limite <- tol_abs + tol_rel * abs(antigo)
    all(abs(novo - antigo) <= limite)
}


#----Extração dos dados----

# Transformações necessárias nos dados:
peru <- read_xlsx("~/Mestrado/Dissertacao/Dados/Inst_Peru.xlsx") %>%
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
    filter(mag >= 4.5)

# Catálogo para pacote ETAS:
peru.quakes <- catalog(peru)

#----Thinning----

# Parâmetros para estimação:
theta_peru <- c("A" = 0.5, "alpha" = 1.0, "c" = 1.0, "p" = 1.15,
                "D" = 0.01, "q" = 1.5, "gamma" = 1.0)

thetas <- list()
thetas[[1]] <- theta_peru

# Número total de eventos:
n_events <- nrow(peru.quakes$revents)

# Definir as sequências do grid explicitamente
x_grid <- seq(range(peru.quakes$region.poly$long)[1], range(peru.quakes$region.poly$long)[2], length.out = 128)
y_grid <- seq(range(peru.quakes$region.poly$lat)[1], range(peru.quakes$region.poly$lat)[2], length.out = 128)

# Raster inicial:
r <- matrix(1, nrow = 128, ncol = 128)
rel_clust <- matrix(0, nrow = 128, ncol = 128)

# Normaliza o raster para a integral ser 1:
# r <- r/(sum(r))

# Na lista de rasters, guardaremos as matrizes:
rasters <- list()
rast_rel_clust <- list()
rasters[[1]] <- r

# Data.frame com os eventos para funções:
df_temp <- data.frame(x = peru.quakes$longlat.coord$long,
                      y = peru.quakes$longlat.coord$lat,
                      mag = peru.quakes$revents[, "mm"] + 4.5,
                      time = peru.quakes$revents[, "tt"])

# Distâncias para kernel adaptativo:
d_j_peru <- d_j(peru.quakes$longlat.coord[, c("long", "lat")])

#----Estimação dos parents/offspring e de background rate----

# Definir o número de núcleos:
num_cores <- 15

t1 <- Sys.time()

passo_e <- calcular_probabilidades(df_temp, theta_peru, 4.5, r, x_grid, y_grid, num_cores)
probabilidades <- passo_e$probabilidades
prob_total <- passo_e$prob_total

# # Paralelização do cálculo de probabilidades:
# probabilidades <- mclapply(1:n_events, function(j) {
#     # Primeiro evento tem probabilidade 0 de ser descendente:
#     if (j == 1) return(numeric(0))
#     # Vetor de índices dos antecessores:
#     i_indices <- 1:(j - 1)
#     # Denominador comum para cada evento j:
#     denom <- lambda_total(
#         df_temp$time[j], df_temp$x[j], df_temp$y[j], df_temp$mag[j],
#         4.5, df_temp, theta_peru, r, x_grid, y_grid
#     )
#     # Computação do numerador para todos os antecessores:
#     numeradores <- kappa_func(df_temp$mag[i_indices], 4.5, theta_peru["A"], theta_peru["alpha"]) *
#         g_func(df_temp$time[j] - df_temp$time[i_indices], theta_peru["c"], theta_peru["p"]) *
#         f_spatial(
#         (df_temp$x[j] - df_temp$x[i_indices])^2 + (df_temp$y[j] - df_temp$y[i_indices])^2,
#         df_temp$mag[i_indices], 4.5,
#         theta_peru["D"], theta_peru["q"], theta_peru["gamma"]
#         )
#     # Retorna o vetor p_{i,j} para o evento j:
#     return(numeradores / denom)
# }, mc.cores = num_cores)
#
# # Probabilidades totais:
# prob_total <- probabilidades %>%
#     lapply(sum) %>%
#     do.call(c, .)

# Estima o background rate:
r <- estimate_background_kernel(x_grid,
                                y_grid,
                                data.frame(x = peru$long,
                                           y = peru$lat,
                                           rho = prob_total),
                                d_j_peru,
                                T_total = max(peru.quakes$rtperiod))

rel_clust <- estimate_relative_clustering(x_grid,
                                          y_grid,
                                          data.frame(x = peru$long,
                                                     y = peru$lat,
                                                     rho = prob_total),
                                          d_j_peru)$C

hessianas <- list()

# Otimiza os parâmetros theta:
otimizacao <- optim(theta_peru,
                        log_lik_etas,
                        gr = gradiente_etas_completo,
                        df_eventos = df_temp,
                        matriz_p_ij = probabilidades,
                        prob_total = prob_total,
                        raster_mu = r, x_grid = x_grid, y_grid = y_grid,
                        method = "L-BFGS-B",
                        lower = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-6),
                        upper = c(15, 10, 15, 15, 15, 15, 10),
                        control = list(maxit = 1000),
                        hessian = T)

theta_peru <- otimizacao$par
thetas[[2]] <- theta_peru

hessianas[[1]] <- otimizacao$hessian

rast_rel_clust[[1]] <- rel_clust

#----Itera até a convergência----

rasters[[2]] <- r

i <- 2

# Trava de segurança:
max_iter_em <- 100

# Tolerâncias de convergência:
tol_rel_raster <- 1e-3
tol_abs_raster <- 1e-6
tol_rel_theta  <- 1e-3
tol_abs_theta  <- 1e-6

# Considerando tolerância de 10e-3:
while (!convergiu_em(rasters[[i]], rasters[[i - 1]], tol_rel_raster, tol_abs_raster) |
       !convergiu_em(thetas[[i]], thetas[[i - 1]], tol_rel_theta, tol_abs_theta)){
    if (i >= max_iter_em) {
        warning(sprintf("EM (ajuste real) atingiu max_iter_em = %d sem convergir dentro da tolerância -- interrompendo.", max_iter_em))
        break
    }
    # Calcula pares de probabilidades (passo E):
    passo_e <- calcular_probabilidades(df_temp, theta_peru, 4.5, r, x_grid, y_grid, num_cores)
    probabilidades <- passo_e$probabilidades
    prob_total <- passo_e$prob_total
    # Estima parâmetros theta:
    otimizacao <- optim(theta_peru,
                        log_lik_etas,
                        gr = gradiente_etas_completo,
                        df_eventos = df_temp,
                        matriz_p_ij = probabilidades,
                        prob_total = prob_total,
                        raster_mu = r, x_grid = x_grid, y_grid = y_grid,
                        method = "L-BFGS-B",
                        lower = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-6),
                        upper = c(15, 10, 15, 15, 15, 15, 10),
                        control = list(maxit = 1000),
                        hessian = T)
    if (otimizacao$convergence != 0) {
        warning(sprintf("optim() não convergiu na iteração %d do ajuste real (código = %d): %s",
                        i, otimizacao$convergence, otimizacao$message))
    }
    theta_peru <- otimizacao$par
    thetas[[i + 1]] <- theta_peru
    hessianas[[i]] <- otimizacao$hessian
    # Estima background rate:
    r <- estimate_background_kernel(x_grid,
                                    y_grid,
                                    data.frame(x = peru$long,
                                               y = peru$lat,
                                               rho = prob_total),
                                    d_j_peru,
                                    T_total = max(peru.quakes$rtperiod))
    rasters[[i + 1]] <- r
    # Estima clustering relativo:
    rel_clust <- estimate_relative_clustering(x_grid,
                                              y_grid,
                                              data.frame(x = peru$long,
                                                         y = peru$lat,
                                                         rho = prob_total),
                                              d_j_peru)$C
    rast_rel_clust[[i]] <- rel_clust
    # Mostra progresso do loop:
    print(i + 1)
    # Atualiza índice:
    i <- i + 1
}

#stopCluster(cl)

t2 <- Sys.time()

print(paste("Tempo total de execução:", t2 - t1))

#----Salvando os Resultados----

save(
    rasters,
    rast_rel_clust,
    thetas,
    probabilidades,
    prob_total,
    hessianas,
    t1, t2,
    file = "resultados_thinning_final_pareto.RData",
    compress = TRUE
)

print("Execução finalizada com sucesso. Arquivos salvos.")
