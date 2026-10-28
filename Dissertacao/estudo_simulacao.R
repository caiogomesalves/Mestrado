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

# Função principal de ajuste com optimParallel:
ajustar_etas_proprio_parallel <- function(
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
                                          cluster = NULL
                                          ) {
    if (is.null(cluster)) {
        return("Especifique o cluster a ser usado na paralelização.")
    }
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
    otim <- optimParallel::optimParallel(theta_atual, log_lik_etas, gr = gradiente_etas_completo,
                                         df_eventos = df_eventos, matriz_p_ij = passo_e$probabilidades,
                                         prob_total = passo_e$prob_total, raster_mu = r_atual,
                                         x_grid = x_grid, y_grid = y_grid, M0 = M0,
                                         method = "L-BFGS-B", lower = lower, upper = upper,
                                         control = list(maxit = 1000), hessian = TRUE,
                                         parallel = list(cl = cluster, forward = F, loginfo = F))
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
        otim <- optimParallel::optimParallel(theta_atual, log_lik_etas, gr = gradiente_etas_completo,
                                             df_eventos = df_eventos, matriz_p_ij = passo_e$probabilidades,
                                             prob_total = passo_e$prob_total, raster_mu = r_atual,
                                             x_grid = x_grid, y_grid = y_grid, M0 = M0,
                                             method = "L-BFGS-B", lower = lower, upper = upper,
                                             control = list(maxit = 1000), hessian = TRUE,
                                             parallel = list(cl = cluster, forward = F, loginfo = F))
        conv_codes  <- c(conv_codes, otim$convergence)
        theta_atual <- otim$par
        hessiana_final <- otim$hessian
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

#----Criar cluster para paralelização----
cl <- makeCluster(15)
setDefaultCluster(cl = cl)
clusterExport(cl, c("d_j", "kappa_func", "g_func", "f_spatial",
                    "func_nu", "lambda_total", "calcular_probabilidades",
                    "log_lik_etas", "gradiente_etas_completo",
                    "estimate_background_kernel", "estimate_relative_clustering",
                    "convergiu_em", "extrair_eventos_catalog"))

#----Ajuste----

region_iran <- list(lat = c(26, 25, 29, 38, 35),
                    long = c(52, 59, 58, 45, 43))

# Catálogo:
iran.cat <- catalog(iran.quakes,
                    study.start = "1991/01/01", study.end = "2011/01/01",
                    region.poly = region_iran, mag.threshold = 4.5)

# Chute inicial parâmetros:
theta_iran <- c("A" = 0.23, "alpha" = 2.8, "c" = 0.022,
                "p" = 1.12, "D" = 0.012, "q" = 2.3, "gamma" = 0.03)

# Definir as sequências do grid explicitamente:
x_grid_iran <- seq(min(iran.cat$longlat.coord$long), max(iran.cat$longlat.coord$long), length.out = 128)
y_grid_iran <- seq(min(iran.cat$longlat.coord$lat), max(iran.cat$longlat.coord$lat), length.out = 128)

t1_teste <- Sys.time()
# Ajuste:
aj_iran <- ajustar_etas_proprio(
    catalogo = iran.cat,
    theta0 = theta_iran,
    M0 = 4.5,
    x_grid = x_grid_iran,
    y_grid = y_grid_iran,
    num_cores = 15,
    guardar_historico = T
)
t2_teste <- Sys.time()

# #====Simulações Irã (bootstrap paramétrico)====
# # Cada réplica: (1) simula um catálogo a partir do modelo final ajustado ao
# # catálogo real (thetas_iran / rasters_iran), (2) reajusta o modelo do zero
# # nesse catálogo simulado pelo mesmo procedimento EM usado no ajuste real.
# # A distribuição das 100 estimativas de theta resultantes é a distribuição
# # bootstrap usada para estimar o erro padrão dos parâmetros.
#
# theta_final_iran   <- thetas_iran[[length(thetas_iran)]]
# raster_final_iran  <- rasters_iran[[length(rasters_iran)]]
#
# # Checagem de criticidade uma única vez: theta_final_iran é o mesmo em todas
# # as réplicas (não muda dentro do loop), então não há motivo para checar (ou
# # silenciar a checagem) 100 vezes.
# verificar_criticidade(theta_final_iran, b = 2.0)
#
# # Grid espacial é fixo (definido pela região de estudo), calculado uma única
# # vez fora do loop em vez de recalculado (de forma idêntica) a cada réplica:
# x_grid_iran <- seq(min(iran.cat$longlat.coord$long), max(iran.cat$longlat.coord$long), length.out = 128)
# y_grid_iran <- seq(min(iran.cat$longlat.coord$lat), max(iran.cat$longlat.coord$lat), length.out = 128)
#
# num_cores    <- 15
# max_iter_em  <- 100  # trava de seguranca por réplica (evita while infinito)
# n_bootstrap  <- 100
# checkpoint_a_cada <- 10  # salva um checkpoint a cada N réplicas concluídas
#
# thetas_ajustados_iran <- vector("list", n_bootstrap)
# backgrounds_iran      <- vector("list", n_bootstrap)
# convergencia_iran     <- vector("list", n_bootstrap)  # diagnostico por replica
# falhas_iran           <- vector("list", n_bootstrap)  # guarda o erro, se houver
#
# for (k in 1:n_bootstrap) {
#     resultado_k <- tryCatch({
#         # Semente distinta e reprodutível por réplica:
#         semente_k <- 10000 + k
#         cat_sim_iran <- simular_catalogo_etas(
#             theta   = theta_final_iran,
#             M0      = 4.5,
#             b       = 2.0,
#             Mmax    = 10,
#             T_max   = 3650,
#             mu_type = "raster",
#             mu_raster = raster_final_iran,
#             x_grid = x_grid_iran,
#             y_grid = y_grid_iran,
#             checar_criticidade = FALSE,
#             semente = semente_k
#         )
#         cat_temp_iran <- catalog(
#             cat_sim_iran %>%
#             mutate(t = ddays(time),
#                    data = POSIXct(1) + t,
#                    long = x, lat = y,
#                    date = date(data),
#                    time = format(data, format = "%H:%M:%S")),
#             region.poly = iran.cat$region.poly
#         )
#         theta_iran_simu <- c("A" = 0.23, "alpha" = 2.8, "c" = 0.022,
#                              "p" = 1.12, "D" = 0.012, "q" = 2.3, "gamma" = 0.03)
#         thetas_iran_temp <- list()
#         thetas_iran_temp[[1]] <- theta_iran_simu
#         # Número total de eventos:
#         n_events_iran_simu <- nrow(cat_temp_iran$revents)
#         # Raster inicial:
#         r_iran_temp <- matrix(1, nrow = 128, ncol = 128)
#         rasters_iran_temp <- list()
#         rasters_iran_temp[[1]] <- r_iran_temp
#         # Data.frame com os eventos para funções:
#         df_temp_iran_temp <- data.frame(x = cat_temp_iran$longlat.coord$long,
#                                         y = cat_temp_iran$longlat.coord$lat,
#                                         mag = cat_temp_iran$revents[, "mm"] + 4.5,
#                                         time = cat_temp_iran$revents[, "tt"])
#         df_temp_iran_temp$time <- df_temp_iran_temp$time - min(df_temp_iran_temp$time)
#         # Distâncias para kernel adaptativo:
#         d_j_iran_temp <- d_j(cat_temp_iran$longlat.coord[, c("long", "lat")])
#         # Passo E inicial (via funcao auxiliar, sem redundancia de calculo):
#         passo_e_temp <- calcular_probabilidades(df_temp_iran_temp, theta_iran_simu, 4.5, r_iran_temp, x_grid_iran, y_grid_iran, num_cores)
#         probabilidades_iran_temp <- passo_e_temp$probabilidades
#         prob_total_iran_temp <- passo_e_temp$prob_total
#         # Estima o background rate:
#         r_iran_temp <- estimate_background_kernel(x_grid_iran,
#                                                   y_grid_iran,
#                                                   data.frame(x = cat_temp_iran$longlat.coord$long,
#                                                              y = cat_temp_iran$longlat.coord$lat,
#                                                              rho = prob_total_iran_temp),
#                                                   d_j_iran_temp,
#                                                   T_total = max(df_temp_iran_temp$time))
#         otimizacao_iran_temp <- optim(theta_iran_simu,
#                                       log_lik_etas,
#                                       gr = gradiente_etas_completo,
#                                       df_eventos = df_temp_iran_temp,
#                                       matriz_p_ij = probabilidades_iran_temp,
#                                       prob_total = prob_total_iran_temp,
#                                       raster_mu = r_iran_temp, x_grid = x_grid_iran, y_grid = y_grid_iran,
#                                       method = "L-BFGS-B",
#                                       lower = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-10),
#                                       upper = c(15, 10, 15, 15, 15, 15, 10),
#                                       control = list(maxit = 1000),
#                                       hessian = T)
#         conv_codes <- otimizacao_iran_temp$convergence
#         theta_iran_simu <- otimizacao_iran_temp$par
#         thetas_iran_temp[[2]] <- theta_iran_simu
#         rasters_iran_temp[[2]] <- r_iran_temp
#         i_em <- 2
#         # Considerando tolerância de 10e-6, com trava de iteracoes maximas:
#         while (!convergiu_em(rasters_iran_temp[[i_em]], rasters_iran_temp[[i_em - 1]], tol_rel_raster, tol_abs_raster) |
#                !convergiu_em(thetas_iran_temp[[i_em]], thetas_iran_temp[[i_em - 1]], tol_rel_theta, tol_abs_theta)){
#                    if (i_em >= max_iter_em) {
#                        warning(sprintf("Replica %d: EM atingiu max_iter_em = %d sem convergir -- interrompendo esta replica.", k, max_iter_em))
#                        break
#                    }
#                    passo_e_temp <- calcular_probabilidades(df_temp_iran_temp, theta_iran_simu, 4.5, r_iran_temp, x_grid_iran, y_grid_iran, num_cores)
#                    probabilidades_iran_temp <- passo_e_temp$probabilidades
#                    prob_total_iran_temp <- passo_e_temp$prob_total
#                    otimizacao_iran_temp <- optim(theta_iran_simu,
#                                                  log_lik_etas,
#                                                  gr = gradiente_etas_completo,
#                                                  df_eventos = df_temp_iran_temp,
#                                                  matriz_p_ij = probabilidades_iran_temp,
#                                                  prob_total = prob_total_iran_temp,
#                                                  raster_mu = r_iran_temp, x_grid = x_grid_iran, y_grid = y_grid_iran,
#                                                  method = "L-BFGS-B",
#                                                  lower = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-10),
#                                                  upper = c(15, 10, 15, 15, 15, 15, 10),
#                                                  control = list(maxit = 1000),
#                                                  hessian = T)
#                    conv_codes <- c(conv_codes, otimizacao_iran_temp$convergence)
#                    theta_iran_simu <- otimizacao_iran_temp$par
#                    thetas_iran_temp[[i_em + 1]] <- theta_iran_simu
#                    r_iran_temp <- estimate_background_kernel(x_grid_iran,
#                                                              y_grid_iran,
#                                                              data.frame(x = cat_temp_iran$longlat.coord$long,
#                                                                         y = cat_temp_iran$longlat.coord$lat,
#                                                                         rho = prob_total_iran_temp),
#                                                              d_j_iran_temp,
#                                                              T_total = max(df_temp_iran_temp$time))
#                    rasters_iran_temp[[i_em + 1]] <- r_iran_temp
#                    cat(sprintf("[Replica %d] iteracao EM %d\n", k, i_em + 1))
#                    i_em <- i_em + 1
#                }
#                list(
#                    theta = thetas_iran_temp[[length(thetas_iran_temp)]],
#                    raster = rasters_iran_temp[[length(rasters_iran_temp)]],
#                    n_iter_em = i_em,
#                    convergence_codes = conv_codes,
#                    convergiu = all(conv_codes == 0) && i_em < max_iter_em
#                )
#     }, error = function(e) {
#         message(sprintf("Replica %d falhou: %s", k, conditionMessage(e)))
#         list(erro = conditionMessage(e))
#     })
#     if (!is.null(resultado_k$erro)) {
#         # Réplica falhou: registra o erro e segue para a próxima, sem
#         # derrubar o job inteiro.
#         falhas_iran[[k]] <- resultado_k$erro
#         thetas_ajustados_iran[[k]] <- NA
#         backgrounds_iran[[k]] <- NA
#         convergencia_iran[[k]] <- NA
#     } else {
#         thetas_ajustados_iran[[k]] <- resultado_k$theta
#         backgrounds_iran[[k]] <- resultado_k$raster
#         convergencia_iran[[k]] <- list(n_iter_em = resultado_k$n_iter_em,
#                                        convergence_codes = resultado_k$convergence_codes,
#                                        convergiu = resultado_k$convergiu)
#     }
#     cat(sprintf("==== Réplica %d/%d concluída ====\n", k, n_bootstrap))
#     # Checkpoint periódico: garante que o trabalho já feito não seja perdido
#     # se o job for interrompido (limite de tempo, falha de nó, etc.) antes
#     # de completar as n_bootstrap réplicas.
#     if (k %% checkpoint_a_cada == 0 || k == n_bootstrap) {
#         save(thetas_ajustados_iran, backgrounds_iran, convergencia_iran, falhas_iran, k,
#              file = "checkpoint_bootstrap_iran.RData", compress = TRUE)
#     }
# }
#
# #----Erro bootstrap dos parâmetros----
#
# # Mantém apenas as réplicas que rodaram sem erro e cujo EM convergiu:
# idx_validas <- which(vapply(convergencia_iran, function(x) !identical(x, NA) && isTRUE(x$convergiu), logical(1)))
#
# if (length(idx_validas) < n_bootstrap) {
#     warning(sprintf("%d de %d réplicas bootstrap foram descartadas (erro na simulação/ajuste ou EM sem convergência). Erro estimado com %d réplicas.",
#                     n_bootstrap - length(idx_validas), n_bootstrap, length(idx_validas)))
# }
#
# theta_matrix_iran <- do.call(rbind, thetas_ajustados_iran[idx_validas])
#
# # Estimativa do erro padrão bootstrap (desvio padrão de cada parâmetro
# # através das réplicas) e do viés bootstrap (média das réplicas menos a
# # estimativa obtida no catálogo real):
# bootstrap_se_iran    <- apply(theta_matrix_iran, 2, sd)
# bootstrap_cov_iran   <- cov(theta_matrix_iran)
# bootstrap_bias_iran  <- colMeans(theta_matrix_iran) - theta_final_iran
# bootstrap_ci_iran    <- apply(theta_matrix_iran, 2, quantile, probs = c(0.025, 0.975))
#
# cat("\n---- Erro padrão bootstrap dos parâmetros ETAS (Irã) ----\n")
# print(bootstrap_se_iran)
# cat("\n---- Viés bootstrap (média das réplicas - estimativa no catálogo real) ----\n")
# print(bootstrap_bias_iran)
# cat("\n---- Intervalo de confiança bootstrap (percentil 2.5%-97.5%) ----\n")
# print(bootstrap_ci_iran)

#----Salvando os Resultados----

save(
    aj_iran,
    # thetas_iran,
    # rasters_iran,
    # rast_rel_clust_iran,
    # thetas_ajustados_iran,
    # backgrounds_iran,
    # probabilidades_iran,
    # prob_total_iran,
    # hessianas_iran,
    # convergencia_iran,
    # falhas_iran,
    # theta_matrix_iran,
    # bootstrap_se_iran,
    # bootstrap_cov_iran,
    # bootstrap_bias_iran,
    # bootstrap_ci_iran,
    file = "teste_ajuste_iran.RData",
    compress = TRUE
)

print("Execução finalizada com sucesso. Arquivos salvos.")
