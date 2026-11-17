#====Pacotes====

library(dplyr)
library(readxl)
library(ETAS)
library(spatstat)
library(parallel)
library(lubridate)

#====Funções====
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
    # S_i é o fator de escala espacial dependente da magnitude:
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

# Calcula as probabilidades p_ij (passo E do EM):
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
        if (length(p_ij) == 0 || sum(p_ij) == 0){
            next
        }
        i_indices <- 1:(j - 1)
        # Diferenças nos tempos:
        dt <- df_eventos$time[j] - df_eventos$time[i_indices]
        # Distâncias espaciais:
        dist2 <- (df_eventos$x[j] - df_eventos$x[i_indices])^2 + (df_eventos$y[j] - df_eventos$y[i_indices])^2
        # Diferenças nas magnitudes:
        dm <- df_eventos$mag[i_indices] - M0
        sum_p_ij <- sum_p_ij + sum(p_ij)
        # Derivada de Produtividade (alpha):
        sum_alpha_term <- sum_alpha_term + sum(p_ij * dm)
        # Derivadas Temporais:
        deriv_c <- ((pp - 1) / cc) - (pp / (dt + cc))
        sum_c_term <- sum_c_term + sum(p_ij * deriv_c)
        deriv_p <- (1 / (pp - 1)) + log(cc) - log(dt + cc)
        sum_p_term <- sum_p_term + sum(p_ij * deriv_p)
        # Derivadas Espaciais (Kernel de Pareto):
        S_i <- D * exp(gamma * dm)
        V_ij <- dist2 / S_i
        U_ij <- 1 + V_ij
        # Fator de escala comum das derivadas em D e gamma:
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
    grad_A     <- - (sum_p_ij / A) + int_A
    grad_alpha <- - sum_alpha_term + int_alpha
    grad_c     <- - sum_c_term
    grad_p     <- - sum_p_term
    grad_D     <- - sum_D_term
    grad_q     <- - sum_q_term
    grad_gamma <- - sum_gamma_term
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
# (no estilo atol + rtol*|valor anterior|).
# Isso evita dois problemas do critério puramente absoluto (|novo - antigo| > 1e-3):
#   (a) para o raster de background, cujos valores costumam ser pequenos
#       (~1e-4 a 1e-2), uma diferença absoluta de 1e-3 já parece "convergida"
#       mesmo representando uma mudança relativa enorme -- isso fazia o loop
#       parar cedo demais, deixando o raster mais próximo do chute inicial
#       uniforme do que deveria (aparência de "suavização excessiva");
#   (b) para os parâmetros theta, que têm escalas muito diferentes entre si
#       (c, D, gamma ~ 0.01-0.03 vs. alpha, p ~ 1-3), a mesma tolerância
#       absoluta é frouxa demais para uns e apertada demais para outros.
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
            "beta (%.3f) <= alpha (%.3f): razão de ramificação é infinita/mal definida.",
            beta, alpha),
            "A distribuição de magnitudes tem cauda mais pesada que o decaimento ",
            "da produtividade -- catálogo pode explodir numericamente.")
        return(invisible(NULL))
    }
    n_bar <- A * beta / (beta - alpha)
    msg <- sprintf("Razão de ramificação estimada: n = %.4f", n_bar)
    if (n_bar >= 1) {
        warning(msg, " -- processo SUPERCRITICO (n >= 1). ",
                "O número esperado de eventos é infinito/instável; a simulação ",
                "provavelmente vai bater no limite de gerações ou demorar muito.")
    } else {
        message(msg, " (subcrítico, processo estável)")
    }
    invisible(n_bar)
}

# Funcao principal de simulação:
simular_catalogo_etas <- function(
                                  theta,                          # Vetor nomeado: A, alpha, c, p, D, q, gamma
                                  M0      = 4.5,                  # Magnitude de corte
                                  b       = 1.0,                  # Valor-b de Gutenberg-Richter
                                  Mmax    = NULL,                 # Magnitude máxima (NULL = sem truncamento)
                                  T_max,                          # Janela temporal de simulação [0, T_max]
                                  mu_type = c("constant", "raster"),
                                  mu_const  = NULL,               # Usado se mu_type == "constant" (eventos/dia/grau^2)
                                  mu_raster = NULL,               # Usado se mu_type == "raster" (matriz nx x ny)
                                  x_grid = NULL, y_grid = NULL,   # Necessários se mu_type == "raster"
                                  xlim = NULL, ylim = NULL,       # Janela espacial (usada se mu_type == "constant")
                                  max_geracoes = 75,              # Trava de segurança contra explosão
                                  checar_criticidade = TRUE,
                                  semente = NULL
                                  ) {
    # Controle de seed para reprodutibilidade:
    if (!is.null(semente)){
        set.seed(semente)
    }
    mu_type <- match.arg(mu_type)
    # Checagem de criticidade do processo:
    if (checar_criticidade){
        verificar_criticidade(theta, b)
    }
    # Parâmetros:
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
        if (n == 0){
            return(numeric(0))
        }
        beta <- b * log(10)
        if (is.null(Mmax)) {
            M0 - log(runif(n)) / beta
        } else {
            u <- runif(n)
            Fmax <- 1 - exp(-beta * (Mmax - M0))
            M0 - log(1 - u * Fmax) / beta
        }
    }
    # ---- Localização de background a partir de um raster (grid) ----
    amostrar_local_raster <- function(n) {
        # Caso nenhuma observação retorna data frame vazio:
        if (n == 0){
            return(data.frame(x = numeric(0),
                              y = numeric(0)))
        }
        n_x <- length(x_grid)
        n_y <- length(y_grid)
        # Probabilidade de cada pixel:
        probs <- as.vector(mu_raster)
        probs <- probs / sum(probs)
        idx <- sample.int(length(probs), size = n, replace = TRUE, prob = probs)
        col_idx <- ((idx - 1) %/% n_x) + 1
        row_idx <- ((idx - 1) %% n_x) + 1
        dx <- diff(x_grid)[1]
        dy <- diff(y_grid)[1]
        x <- x_grid[row_idx] + runif(n, -dx / 2, dx / 2)
        y <- y_grid[col_idx] + runif(n, -dy / 2, dy / 2)
        data.frame(x = x, y = y)
    }
    # ---- Deslocamento espacial de réplicas: kernel de Pareto isotrópico ----
    # CDF do raio (obtida integrando f(x,y) em coordenadas polares):
    #   F(r) = 1 - (1 + r^2/S)^{-(q-1)}   =>   r^2 = S * [(1-u)^{-1/(q-1)} - 1]
    simular_offset_espacial <- function(mag) {
        n <- length(mag)
        if (n == 0) {
            return(data.frame(dx = numeric(0),
                              dy = numeric(0)))
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
        if (n == 0){
            return(numeric(0))
        }
        u <- runif(n)
        cc * ((1 - u)^(-1 / (pp - 1)) - 1)
    }
    # ---- Passo 1: Eventos de background ----
    # Caso raster não seja passado:
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
        # Número de observações amostradas do background:
        n_bg <- rpois(1, mu_total * T_max)
        # Localizações dos pontos:
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
    # ---- Passo 2: Processar réplicas geração por geração (fila / BFS) ----
    proximo_id <- n_bg + 1L
    fila <- catalogo
    g_final <- 0L
    if (n_bg > 0) {
        for (g in seq_len(max_geracoes)) {
            g_final <- g
            # Critério de parada (todas as gerações finalizaram):
            if (nrow(fila) == 0){
                break
            }
            # Produtividade:
            kappa_vals <- A * exp(alpha * (fila$mag - M0))
            # Número de filhos produzidos:
            n_filhos   <- rpois(nrow(fila), kappa_vals)
            if (sum(n_filhos) == 0){
                break
            }
            idx_pais <- rep(seq_len(nrow(fila)), times = n_filhos)
            pais <- fila[idx_pais, , drop = FALSE]
            # Magnitudes dos filhos:
            mags_filhos <- simular_magnitudes(nrow(pais))
            # Diferenças temporais dos filhos:
            dt          <- simular_offset_temporal(nrow(pais))
            # Distâncias espaciais dos filhos:
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
            # Descarta filhos fora da janela temporal:
            filhos <- filhos[filhos$time <= T_max, , drop = FALSE]
            catalogo <- rbind(catalogo, filhos)
            fila <- filhos
        }
        if (g_final == max_geracoes && nrow(fila) > 0) {
            warning("Número máximo de gerações atingido -- verifique a razão de ",
                    "ramificação (verificar_criticidade) antes de confiar no catálogo.")
        }
    }
    # Ordena o catálogo temporalmente:
    catalogo <- catalogo[order(catalogo$time), ]
    rownames(catalogo) <- NULL
    # Atributos invisíveis do catálogo:
    attr(catalogo, "theta_verdadeiro") <- theta
    attr(catalogo, "M0") <- M0
    attr(catalogo, "b")  <- b
    catalogo
}

# Função principal de ajuste:
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
                                 guardar_historico = FALSE
                                 ) {
    df_eventos <- extrair_eventos_catalog(catalogo, M0)
    df_eventos$time <- df_eventos$time - min(df_eventos$time)
    if (is.null(d_j_pre)) d_j_pre <- d_j(df_eventos[, c("x", "y")])
    if (is.null(r0))      r0      <- matrix(1, nrow = length(x_grid), ncol = length(y_grid))
    thetas_hist  <- list(theta0)
    rasters_hist <- list(r0)
    theta_atual  <- theta0
    r_atual      <- r0
    conv_codes   <- integer(0)
    # ---- Primeira iteração do EM (sempre executada, sem checagem prévia) ----
    passo_e <- calcular_probabilidades(df_eventos, theta_atual, M0, r_atual, x_grid, y_grid, num_cores)
    r_atual <- estimate_background_kernel(x_grid, y_grid,
                                          data.frame(x = df_eventos$x, y = df_eventos$y, rho = passo_e$prob_total),
                                          d_j_pre, T_total = max(df_eventos$time))
    otim <- optim(theta_atual, log_lik_etas, gr = gradiente_etas_completo,
                  df_eventos = df_eventos, matriz_p_ij = passo_e$probabilidades,
                  prob_total = passo_e$prob_total, raster_mu = r_atual,
                  x_grid = x_grid, y_grid = y_grid, M0 = M0,
                  method = "L-BFGS-B", lower = lower, upper = upper,
                  control = list(maxit = 1000), hessian = TRUE)
    conv_codes  <- c(conv_codes, otim$convergence)
    theta_atual <- otim$par
    hessiana_final <- otim$hessian
    thetas_hist[[2]]  <- theta_atual
    rasters_hist[[2]] <- r_atual
    # ---- Iterações seguintes, até convergir ou atingir max_iter_em ----
    i <- 2
    repeat {
        if (convergiu_em(rasters_hist[[i]], rasters_hist[[i - 1]], tol_rel_raster, tol_abs_raster) &&
            convergiu_em(thetas_hist[[i]], thetas_hist[[i - 1]], tol_rel_theta, tol_abs_theta)) {
            break
        }
        if (i >= max_iter_em) {
            warning(sprintf("ajustar_etas_proprio: max_iter_em = %d atingido sem convergir.", max_iter_em))
            break
        }
        passo_e <- calcular_probabilidades(df_eventos, theta_atual, M0, r_atual, x_grid, y_grid, num_cores)
        r_atual <- estimate_background_kernel(x_grid, y_grid,
                                              data.frame(x = df_eventos$x, y = df_eventos$y, rho = passo_e$prob_total),
                                              d_j_pre, T_total = max(df_eventos$time))
        otim <- optim(theta_atual, log_lik_etas, gr = gradiente_etas_completo,
                      df_eventos = df_eventos, matriz_p_ij = passo_e$probabilidades,
                      prob_total = passo_e$prob_total, raster_mu = r_atual,
                      x_grid = x_grid, y_grid = y_grid, M0 = M0,
                      method = "L-BFGS-B", lower = lower, upper = upper,
                      control = list(maxit = 1000), hessian = TRUE)
        conv_codes  <- c(conv_codes, otim$convergence)
        theta_atual <- otim$par
        hessiana_final <- otim$hessian
        cat("====Iteração", i, "finalizada====")
        cat(sprintf("\n"))
        i <- i + 1
        thetas_hist[[i]]  <- theta_atual
        rasters_hist[[i]] <- r_atual
    }
    list(
        theta              = theta_atual,
        raster             = r_atual,
        hessiana           = hessiana_final,
        n_iter_em          = i,
        convergence_codes  = conv_codes,
        convergiu          = all(conv_codes == 0) && i < max_iter_em,
        thetas_hist        = if (guardar_historico) thetas_hist  else NULL,
        rasters_hist       = if (guardar_historico) rasters_hist else NULL
    )
}

# Converte um objeto catalog (pacote ETAS) ou data.frame:
extrair_eventos_catalog <- function(catalogo_historico, M0) {
    if (inherits(catalogo_historico, "catalog")) {
        idx <- catalogo_historico$revents[, "flag"] %in% c(0, 1)
        df <- data.frame(
            time = catalogo_historico$revents[idx, "tt"],
            x    = catalogo_historico$longlat.coord$long[idx],
            y    = catalogo_historico$longlat.coord$lat[idx],
            mag  = catalogo_historico$revents[idx, "mm"] + M0
        )
    } else if (is.data.frame(catalogo_historico)) {
        obrigatorias <- c("time", "x", "y", "mag")
        faltando <- setdiff(obrigatorias, names(catalogo_historico))
        if (length(faltando) > 0) {
            stop(sprintf("catalogo_historico está sem a(s) coluna(s): %s", paste(faltando, collapse = ", ")))
        }
        df <- catalogo_historico[, obrigatorias, drop = FALSE]
    } else {
        stop("catalogo_historico deve ser um data.frame (colunas time, x, y, mag) ou um objeto 'catalog' do pacote ETAS.")
    }
    if (!"id" %in% names(df)) df$id <- seq_len(nrow(df))
    df <- df[order(df$time), ]
    rownames(df) <- NULL
    df
}

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

#====Extração dos dados====

# Transformações necessárias nos dados:
peru <- read_xlsx("Inst_Peru.xlsx") %>%
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

#====Ajuste====

# Parâmetros para estimação:
theta_peru <- c("A" = 0.5, "alpha" = 1.0, "c" = 1.0, "p" = 1.15,
                "D" = 0.01, "q" = 1.5, "gamma" = 1.0)

# Sequências do grid:
x_grid <- seq(range(peru.quakes$region.poly$long)[1], range(peru.quakes$region.poly$long)[2], length.out = 64)
y_grid <- seq(range(peru.quakes$region.poly$lat)[1], range(peru.quakes$region.poly$lat)[2], length.out = 64)

t1_peru <- Sys.time()

# Ajuste:
aj_peru <- ajustar_etas_proprio(
    catalogo = peru.quakes,
    theta0 = theta_peru,
    M0 = 4.5,
    x_grid = x_grid,
    y_grid = y_grid,
    num_cores = 15,
    guardar_historico = T
)

t2_peru <- Sys.time()

#====Salvando os Resultados====

save(
    aj_peru,
    t1_peru,
    t2_peru,
    file = "ajuste_peru.RData",
    compress = TRUE
)
