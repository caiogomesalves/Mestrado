#----Pacotes----

library(dplyr)
library(ETAS)
library(spatstat)
library(parallel)
library(lubridate)

#----Funções----

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

# Calcula as probabilidades p_ij (passo E do EM) para todos os eventos de uma vez.
# Substitui o padrão anterior (chamar lambda_total() para obter o denominador e,
# em seguida, recalcular manualmente os mesmos numeradores para a soma de
# clustering): aqui o denominador é obtido diretamente como mu_j + sum(numeradores),
# evitando refazer o mesmo somatório duas vezes (uma via dplyr dentro de
# lambda_total, outra vetorizada) a cada eventos e a cada iteração do EM.
# Também evita chamar func_nu() ponto a ponto: os valores de mu são extraídos
# de uma só vez, vetorizados, para todos os eventos.
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

# Critério de convergência do EM: combina tolerância absoluta e relativa
# (no estilo atol + rtol*|valor anterior|, usado em solvers numéricos).
# Isso evita dois problemas do critério puramente absoluto (|novo - antigo| > 1e-3):
#   (a) para o raster de background, cujos valores costumam ser pequenos
#       (~1e-4 a 1e-2), uma diferença absoluta de 1e-3 já parece "convergida"
#       mesmo representando uma mudança relativa enorme -- isso fazia o loop
#       parar cedo demais, deixando o raster mais próximo do chute inicial
#       uniforme do que deveria (aparência de "suavização excessiva");
#   (b) para os parâmetros theta, que têm escalas muito diferentes entre si
#       (c, D, gamma ~ 0.01-0.03 vs. alpha, p ~ 1-3), a mesma tolerância
#       absoluta é frouxa demais para uns e apertada demais para outros.
# tol_abs funciona como piso para valores próximos de zero, evitando que a
# parte relativa (rtol*|antigo|) exploda por divisão por valor pequeno.
convergiu_em <- function(novo, antigo, tol_rel = 1e-3, tol_abs = 1e-6) {
    limite <- tol_abs + tol_rel * abs(antigo)
    all(abs(novo - antigo) <= limite)
}

# Verifica criticidade do processo:
verificar_criticidade <- function(theta, b) {
  A <- unname(theta["A"]); alpha <- unname(theta["alpha"])
  beta <- b * log(10)
  if (beta <= alpha) {
    warning(sprintf(
      "beta (%.3f) <= alpha (%.3f): razao de ramificacao e' infinita/mal definida. ",
      beta, alpha),
      "A distribuicao de magnitudes tem cauda mais pesada que o decaimento ",
      "da produtividade -- catalogo pode explodir numericamente.")
    return(invisible(NULL))
  }
  n_bar <- A * beta / (beta - alpha)
  msg <- sprintf("Razao de ramificacao estimada: n = %.4f", n_bar)
  if (n_bar >= 1) {
    warning(msg, " -- processo SUPERCRITICO (n >= 1). ",
            "O numero esperado de eventos e' infinito/instavel; a simulacao ",
            "provavelmente vai bater no limite de geracoes ou demorar muito.")
  } else {
    message(msg, " (subcritico, processo estavel)")
  }
  invisible(n_bar)
}

# Funcao principal de simulação:
simular_catalogo_etas <- function(
                                  theta,                          # vetor nomeado: A, alpha, c, p, D, q, gamma
                                  M0      = 4.5,                  # magnitude de corte (completude)
                                  b       = 1.0,                  # valor-b de Gutenberg-Richter p/ magnitudes
                                  Mmax    = NULL,                 # Magnitude maxima (NULL = sem truncamento)
                                  T_max,                          # Janela temporal de simulacao [0, T_max]
                                  mu_type = c("constant", "raster"),
                                  mu_const  = NULL,               # Usado se mu_type == "constant" (eventos/dia/grau^2)
                                  mu_raster = NULL,               # Usado se mu_type == "raster" (matriz nx x ny)
                                  x_grid = NULL, y_grid = NULL,   # Necessarios se mu_type == "raster"
                                  xlim = NULL, ylim = NULL,       # Janela espacial (usada se mu_type == "constant")
                                  max_geracoes = 75,              # Trava de seguranca contra explosão
                                  checar_criticidade = TRUE,
                                  semente = NULL
                                  ) {
    if (!is.null(semente)){
        set.seed(semente)
    }
    mu_type <- match.arg(mu_type)
    if (checar_criticidade){
        verificar_criticidade(theta, b)
    }
    A     <- unname(theta["A"])
    alpha <- unname(theta["alpha"])
    cc    <- unname(theta["c"])
    pp    <- unname(theta["p"])
    D     <- unname(theta["D"])
    q     <- unname(theta["q"])
    gamma <- unname(theta["gamma"])
    # Garante que valores sejam próprios para densidades:
    stopifnot(pp > 1, q > 1)
    # ---- Magnitudes: Gutenberg-Richter (truncada se Mmax fornecido) ----
    simular_magnitudes <- function(n) {
        if (n == 0) return(numeric(0))
        beta <- b * log(10)
        if (is.null(Mmax)) {
            M0 - log(runif(n)) / beta
        } else {
            u <- runif(n)
            Fmax <- 1 - exp(-beta * (Mmax - M0))
            M0 - log(1 - u * Fmax) / beta
        }
    }
    # ---- Localizacao de background a partir de um raster (grid) ----
    amostrar_local_raster <- function(n) {
        # Caso nenhuma observação retorna data frame vazio:
        if (n == 0){
            return(data.frame(x = numeric(0), y = numeric(0)))
        }
        n_x <- length(x_grid)
        n_y <- length(y_grid)
        probs <- as.vector(mu_raster)
        probs <- probs / sum(probs)
        idx <- sample.int(length(probs), size = n, replace = TRUE, prob = probs)
        col_idx <- ((idx - 1) %/% n_x) + 1
        row_idx <- ((idx - 1) %% n_x) + 1
        dx <- diff(x_grid)[1]; dy <- diff(y_grid)[1]
        x <- x_grid[row_idx] + runif(n, -dx / 2, dx / 2)
        y <- y_grid[col_idx] + runif(n, -dy / 2, dy / 2)
        data.frame(x = x, y = y)
    }
    # ---- Deslocamento espacial de replicas: kernel de Pareto isotropico ----
    # CDF do raio (obtida integrando f(x,y) em coordenadas polares):
    #   F(r) = 1 - (1 + r^2/S)^{-(q-1)}   =>   r^2 = S * [(1-u)^{-1/(q-1)} - 1]
    simular_offset_espacial <- function(mag) {
        n <- length(mag)
        if (n == 0) {
            return(data.frame(dx = numeric(0), dy = numeric(0)))
        }
        S_i <- D * exp(gamma * (mag - M0))
        u   <- runif(n)
        r2  <- S_i * ((1 - u)^(-1 / (q - 1)) - 1)
        r   <- sqrt(r2)
        ang <- runif(n, 0, 2 * pi)
        data.frame(dx = r * cos(ang), dy = r * sin(ang))
    }
    # ---- Deslocamento temporal: Omori-Utsu (CDF inversa) ----
    # F(t) = 1 - (c/(t+c))^{p-1}   =>   t = c * [(1-u)^{-1/(p-1)} - 1]
    simular_offset_temporal <- function(n) {
        if (n == 0) return(numeric(0))
        u <- runif(n)
        cc * ((1 - u)^(-1 / (pp - 1)) - 1)
    }
    # ---- 1. Eventos de background (geracao 0) ----
    if (mu_type == "constant") {
        stopifnot(!is.null(mu_const), !is.null(xlim), !is.null(ylim))
        area <- diff(xlim) * diff(ylim)
        n_bg <- rpois(1, mu_const * area * T_max)
        x0 <- runif(n_bg, xlim[1], xlim[2])
        y0 <- runif(n_bg, ylim[1], ylim[2])
    } else {
        stopifnot(!is.null(mu_raster), !is.null(x_grid), !is.null(y_grid))
        dxg <- diff(x_grid)[1]
        dyg <- diff(y_grid)[1]
        # Integral aproximada de mu(x,y):
        mu_total <- sum(mu_raster) * dxg * dyg
        n_bg <- rpois(1, mu_total * T_max)
        loc <- amostrar_local_raster(n_bg)
        x0 <- loc$x
        y0 <- loc$y
    }
    t0 <- runif(n_bg, 0, T_max)
    m0 <- simular_magnitudes(n_bg)
    catalogo <- data.frame(
        id = seq_len(n_bg), parent = 0L, geracao = 0L,
        time = t0, x = x0, y = y0, mag = m0
    )
    # ---- 2. Processar réplicas geração por geração (fila / BFS) ----
    proximo_id <- n_bg + 1L
    fila <- catalogo
    g_final <- 0L
    if (n_bg > 0) {
        for (g in seq_len(max_geracoes)) {
            g_final <- g
            if (nrow(fila) == 0){
                break
            }
            kappa_vals <- A * exp(alpha * (fila$mag - M0))
            n_filhos   <- rpois(nrow(fila), kappa_vals)
            if (sum(n_filhos) == 0){
                break
            }
            idx_pais <- rep(seq_len(nrow(fila)), times = n_filhos)
            pais <- fila[idx_pais, , drop = FALSE]
            mags_filhos <- simular_magnitudes(nrow(pais))
            dt          <- simular_offset_temporal(nrow(pais))
            off_esp     <- simular_offset_espacial(pais$mag)
            filhos <- data.frame(
                id      = proximo_id + seq_len(nrow(pais)) - 1L,
                parent  = pais$id,
                geracao = g,
                time    = pais$time + dt,
                x       = pais$x + off_esp$dx,
                y       = pais$y + off_esp$dy,
                mag     = mags_filhos
            )
            proximo_id <- proximo_id + nrow(filhos)
            # Descarta filhos fora da janela temporal
            filhos <- filhos[filhos$time <= T_max, , drop = FALSE]
            catalogo <- rbind(catalogo, filhos)
            fila <- filhos
        }
        if (g_final == max_geracoes && nrow(fila) > 0) {
            warning("Número máximo de gerações atingido -- verifique a razão de ",
                    "ramificação (verificar_criticidade) antes de confiar no catálogo.")
        }
    }
    catalogo <- catalogo[order(catalogo$time), ]
    rownames(catalogo) <- NULL
    attr(catalogo, "theta_verdadeiro") <- theta
    attr(catalogo, "M0") <- M0
    attr(catalogo, "b")  <- b
    catalogo
}

#----Ajuste----

sumatra_quakes <- read.csv("~/Mestrado/Dissertacao/Dados/sumatra_quakes.csv")

sumatra_cat <- sumatra_quakes %>%
    arrange(desc(time)) %>%
    select(date = time,
           time = time,
           lat = latitude, long = longitude,
           prof = depth,
           mag) %>%
    mutate(
        date = as.Date(date),
        time = (as.POSIXct(time, format="%Y-%m-%dT%H:%M:%OSZ", tz = "UTC") %>%
                format("%H:%M:%S")),
        lat = as.numeric(lat),
        long = as.numeric(long),
        prof = as.numeric(prof),
        mag = as.numeric(mag)
    ) %>%
    catalog()

# Chute inicial parâmetros:
theta_sum <- c("A" = 0.01, "alpha" = 1, "c" = 0.01,
               "p" = 1.3, "D" = 0.01, "q" = 2, "gamma" = 1)

thetas_sum <- list()
thetas_sum[[1]] <- theta_sum

# Número total de eventos:
n_events_sum <- nrow(sumatra_cat$revents)

# Definir as sequências do grid explicitamente
x_grid_sum <- seq(min(sumatra_cat$longlat.coord$long), max(sumatra_cat$longlat.coord$long), length.out = 128)
y_grid_sum <- seq(min(sumatra_cat$longlat.coord$lat), max(sumatra_cat$longlat.coord$lat), length.out = 128)

# Raster inicial:
r_sum <- matrix(1, nrow = 128, ncol = 128)
rel_clust_sum <- matrix(0, nrow = 128, ncol = 128)

# Normaliza o raster para a integral ser 1:
# r_iran <- r_iran/(sum(r_iran))

# Na lista de rasters, guardaremos as matrizes:
rasters_sum <- list()
rast_rel_clust_sum <- list()
rasters_sum[[1]] <- r_sum

hessianas_sum <- list()

# Data.frame com os eventos para funções:
df_temp_sum <- data.frame(x = sumatra_cat$longlat.coord$long,
                          y = sumatra_cat$longlat.coord$lat,
                          mag = sumatra_cat$revents[, "mm"] + 4,
                          time = sumatra_cat$revents[, "tt"])

# Distâncias para kernel adaptativo:
d_j_sum <- d_j(sumatra_cat$longlat.coord[, c("long", "lat")])

num_cores <- 15

t1_sum <- Sys.time()

# Cálculo (paralelo) das probabilidades p_ij e das probabilidades totais:
passo_e_sum <- calcular_probabilidades(df_temp_sum, theta_sum, 4, r_sum, x_grid_sum, y_grid_sum, num_cores)
probabilidades_sum <- passo_e_sum$probabilidades
prob_total_sum <- passo_e_sum$prob_total

# Estima o background rate:
r_sum <- estimate_background_kernel(x_grid_sum,
                                    y_grid_sum,
                                    data.frame(x = sumatra_cat$longlat.coord$long,
                                               y = sumatra_cat$longlat.coord$lat,
                                               rho = prob_total_sum),
                                    d_j_sum,
                                    T_total = max(df_temp_sum$time))

rel_clust_sum <- estimate_relative_clustering(x_grid_sum,
                                              y_grid_sum,
                                              data.frame(x = sumatra_cat$longlat.coord$long,
                                                         y = sumatra_cat$longlat.coord$lat,
                                                         rho = prob_total_sum),
                                              d_j_sum)$C

# Otimiza os parâmetros theta:
otimizacao_sum <- optim(theta_sum,
                        log_lik_etas,
                        gr = gradiente_etas_completo,
                        df_eventos = df_temp_sum,
                        matriz_p_ij = probabilidades_sum,
                        prob_total = prob_total_sum,
                        raster_mu = r_sum, M0 = 4,
                        x_grid = x_grid_sum, y_grid = y_grid_sum,
                        method = "L-BFGS-B",
                        lower = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-10),
                        upper = c(15, 10, 15, 15, 15, 15, 10),
                        control = list(maxit = 1000),
                        hessian = T)

theta_sum <- otimizacao_sum$par
thetas_sum[[2]] <- theta_sum

hessianas_sum[[1]] <- otimizacao_sum$hessian

rast_rel_clust_sum[[1]] <- rel_clust_sum

rasters_sum[[2]] <- r_sum

i <- 2

# Trava de seguranca: numero maximo de iteracoes do EM (evita loop infinito
# caso a convergencia oscile e nunca fique abaixo da tolerancia):
max_iter_em <- 100

# Tolerâncias de convergência (ajustáveis). tol_rel é a mudança relativa
# máxima tolerada entre iterações (ex.: 1e-3 = 0.1%); tol_abs é o piso
# absoluto usado para valores próximos de zero, para evitar que a parte
# relativa fique instável nesse regime.
tol_rel_raster <- 1e-3
tol_abs_raster <- 1e-6
tol_rel_theta  <- 1e-3
tol_abs_theta  <- 1e-6

# Considerando tolerância de 10e-6:
while (!convergiu_em(rasters_sum[[i]], rasters_sum[[i - 1]], tol_rel_raster, tol_abs_raster) |
       !convergiu_em(thetas_sum[[i]], thetas_sum[[i - 1]], tol_rel_theta, tol_abs_theta)){
    if (i >= max_iter_em) {
        warning(sprintf("EM (ajuste real) atingiu max_iter_em = %d sem convergir dentro da tolerancia -- interrompendo.", max_iter_em))
        break
    }
    # Calcula pares de probabilidades (passo E):
    passo_e_sum <- calcular_probabilidades(df_temp_sum, theta_sum, 4, r_sum, x_grid_sum, y_grid_sum, num_cores)
    probabilidades_sum <- passo_e_sum$probabilidades
    prob_total_sum <- passo_e_sum$prob_total
    # Estima parâmetros theta:
    otimizacao_sum <- optim(theta_sum,
                            log_lik_etas,
                            gr = gradiente_etas_completo,
                            df_eventos = df_temp_sum,
                            matriz_p_ij = probabilidades_sum,
                            prob_total = prob_total_sum,
                            raster_mu = r_sum, M0 = 4,
                            x_grid = x_grid_sum, y_grid = y_grid_sum,
                            method = "L-BFGS-B",
                            lower = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-10),
                            upper = c(15, 10, 15, 15, 15, 15, 10),
                            control = list(maxit = 1000),
                            hessian = T)
    if (otimizacao_sum$convergence != 0) {
        warning(sprintf("optim() nao convergiu na iteracao %d do ajuste real (codigo = %d): %s",
                        i, otimizacao_sum$convergence, otimizacao_sum$message))
    }
    theta_sum <- otimizacao_sum$par
    thetas_sum[[i + 1]] <- theta_sum
    hessianas_sum[[i]] <- otimizacao_sum$hessian
    # Estima background rate:
    r_sum <- estimate_background_kernel(x_grid_sum,
                                        y_grid_sum,
                                        data.frame(x = sumatra_cat$longlat.coord$long,
                                                   y = sumatra_cat$longlat.coord$lat,
                                                   rho = prob_total_sum),
                                        d_j_sum,
                                        T_total = max(df_temp_sum$time))
    rasters_sum[[i + 1]] <- r_sum
    # Estima clustering relativo:
    rel_clust_sum <- estimate_relative_clustering(x_grid_sum,
                                                  y_grid_sum,
                                                  data.frame(x = sumatra_cat$longlat.coord$long,
                                                             y = sumatra_cat$longlat.coord$lat,
                                                             rho = prob_total_sum),
                                                  d_j_sum)$C
    rast_rel_clust_sum[[i]] <- rel_clust_sum
    # Mostra progresso do loop:
    print(i + 1)
    # Atualiza índice:
    i <- i + 1
}

t2_sum <- Sys.time()

print(paste("Tempo total de execução:", t2_sum - t1_sum))

#----Salvando os Resultados----

# save(
#     rasters_sum,
#     rast_rel_clust_sum,
#     thetas_sum,
#     probabilidades_sum,
#     prob_total_sum,
#     hessianas_sum,
#     t1_sum, t2_sum,
#     file = "resultados_ajuste_sumatra.RData",
#     compress = TRUE
# )

print("Execução finalizada com sucesso. Arquivos salvos.")

#----Resíduos----

load("~/resultados_ajuste_sumatra.RData")

# Taxa global de background (Eventos/dia em toda a região):
mu_global_sum <- sum(pmax(0, 1 - prob_total_sum)) / max(df_temp_sum$time)

# Tau_j para os tempos dos eventos:
tau_j_sum <- calcular_tau(df_temp_sum$time, df_temp_sum, thetas_sum[[length(thetas_sum)]], mu_global_sum, M0 = 4)

# Calcular U_j:
delta_tau_sum <- diff(c(0, tau_j_sum))
U_j_sum <- 1 - exp(-delta_tau_sum)

#----Gráficos----
par(mfrow = c(1, 2))

# Gráfico A: tau_j vs j
plot(1:length(tau_j_sum), tau_j_sum, type = "l", col = "blue", lwd = 2,
     xlab = "Número do Evento (j)", ylab = expression(tau[j]),
     main = "Tempos Transformados")
abline(a = 0, b = 1, col = "red", lty = 2)

# Gráfico B: Q-Q Plot de U_j
plot(sort(U_j_sum), qunif(ppoints(length(U_j_sum))), pch = 20, col = "darkgreen",
     xlab = "Empírico", ylab = "Teórico (Uniforme 0-1)",
     main = "Q-Q Plot de U_j")
abline(a = 0, b = 1, col = "red", lty = 2)

ks.test(U_j_sum, punif)

# Número de caixas (bins) de tempo
n_temp_sum <- 100
janelas_tempo_sum <- seq(0, max(df_temp_sum$time), length.out = n_temp_sum + 1)

# Contagem de eventos em cada janela:
obs_temp_sum <- hist(df_temp_sum$time, breaks = janelas_tempo_sum, plot = FALSE)$counts

# Contagens esperadas pela diferença do tau nos limites da janela:
tau_limites_sum <- calcular_tau(janelas_tempo_sum, df_temp_sum, thetas_sum[[length(thetas_sum)]], mu_global_sum, M0 = 4)
esp_temp_sum <- diff(tau_limites_sum)

# Resíduos de Pearson Temporais:
res_temp_pearson_sum <- (obs_temp_sum - esp_temp_sum) / sqrt(esp_temp_sum)

# Resíduos Temporais:
par(mfrow = c(1, 1))
plot(janelas_tempo_sum[-1], res_temp_pearson_sum, type = "h", lwd = 2, col = "purple",
     xlab = "Tempo (dias)", ylab = "Resíduos de Pearson",
     main = "Resíduos Temporais Suavizados")
abline(h = 0, col = "black")
abline(h = c(-2, 2), col = "red", lty = 2)

#====Espaciais====

# Definir os pontos observados:
pontos_sum <- ppp(df_temp_sum$x, df_temp_sum$y,
                  window = owin(x_grid_sum[c(1, length(x_grid_sum))],
                                y_grid_sum[c(1, length(y_grid_sum))]),
                  marks = pmax(0, 1 - prob_total_sum))

# Criar a densidade espacial esperada pelo Background:
cell_area_sum <- diff(x_grid_sum)[1] * diff(y_grid_sum)[1]
matriz_esperada_sum <- (rasters_sum[[length(rasters_sum)]] / (sum(rasters_sum[[length(rasters_sum)]]) * cell_area_sum)) * (mu_global_sum * max(df_temp_sum$time))
# matriz_esperada_iran <- (rasters_iran[[length(rasters_iran)]] / sum(rasters_iran[[length(rasters_iran)]])) * (mu_global_iran * max(df_temp_iran$time))

# Extrair os 7 parâmetros do modelo convergido
A_sum <- thetas_sum[[length(thetas_sum)]]["A"]
alpha_sum <- thetas_sum[[length(thetas_sum)]]["alpha"]
D_sum <- thetas_sum[[length(thetas_sum)]]["D"]
q_par_sum <- thetas_sum[[length(thetas_sum)]]["q"]
gamma_par_sum <- thetas_sum[[length(thetas_sum)]]["gamma"]

M0_sum <- 4 # Ajuste para o limiar de magnitude correto do seu catálogo atual

# Sobrepor a nuvem espacial de Pareto (Lei de Potência) de TODOS os eventos
for(i in 1:nrow(df_temp_sum)) {
    xi_sum <- df_temp_sum$x[i]
    yi_sum <- df_temp_sum$y[i]
    mi_sum <- df_temp_sum$mag[i]
    # Produtividade esperada do evento (mantém-se com alpha)
    k_i_sum <- A_sum * exp(alpha_sum * (mi_sum - M0_sum))
    # Fator de escala espacial S_i (agora independente de alpha, usa gamma)
    S_i_sum <- D_sum * exp(gamma_par_sum * (mi_sum - M0_sum))
    # Matriz de distâncias ao quadrado
    dx2_sum <- (x_grid_sum - xi_sum)^2
    dy2_sum <- (y_grid_sum - yi_sum)^2
    dist2_matrix_sum <- outer(dx2_sum, dy2_sum, "+")
    # Kernel de Pareto Bivariado Isotrópico (Fórmula do pacote ETAS)
    # f(x,y) = ((q-1) / (pi * S_i)) * (1 + dist2 / S_i)^(-q)
    term1_sum <- (q_par_sum - 1) / (pi * S_i_sum)
    term2_sum <- (1 + (dist2_matrix_sum / S_i_sum))^(-q_par_sum)
    kernel_esp_sum <- k_i_sum * term1_sum * term2_sum
    matriz_esperada_sum <- matriz_esperada_sum + kernel_esp_sum
}

densidade_esp_sum <- im(t(matriz_esperada_sum), xcol = x_grid_sum, yrow = y_grid_sum)

# Criar a densidade observada, forçando o alinhamento com xy:
densidade_obs_sum <- density(pontos_sum, sigma = bw.diggle(pontos_sum), xy = densidade_esp_sum,
                             weights = pontos_sum$marks)

# Calcular os resíduos de Pearson:
epsilon <- 1e-10
res_spat_pearson_sum <- eval.im((densidade_obs_sum - densidade_esp_sum) / sqrt(densidade_esp_sum + epsilon))

# Visualização:
plot(res_spat_pearson_sum, main = "Resíduos Espaciais Totais de Pearson", col = topo.colors(256))
