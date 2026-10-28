#====Funções auxiliares====

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
                                 catalogo,                 # data.frame (x,y,mag,time) OU objeto catalog do pacote ETAS
                                 theta0,                   # vetor nomeado: chute inicial (A, alpha, c, p, D, q, gamma)
                                 M0 = 4.5,
                                 x_grid, y_grid,
                                 d_j_pre = NULL,            # distâncias pré-computadas (evita recálculo entre chamadas com o mesmo catálogo)
                                 r0 = NULL,                 # raster inicial de fundo (default: uniforme)
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

#====Resíduos====

# Função para os Tempos Transformados (tau):
calcular_tau <- function(tempos_eval, df_cat, params, mu_global, M0 = 4) {
    A <- params["A"]
    alpha <- params["alpha"]
    cc <- params["c"]
    pp <- params["p"]
    tau_vetor <- numeric(length(tempos_eval))
    for (k in 1:length(tempos_eval)) {
        t_k <- tempos_eval[k]
        # Integral do Background (Linear no tempo):
        tau_bg <- mu_global * t_k
        # Integral do Clustering (Lei de Omori)
        idx <- which(df_cat$time < t_k)
        if (length(idx) > 0) {
            dt <- t_k - df_cat$time[idx]
            kappa_vals <- A * exp(alpha * (df_cat$mag[idx] - M0))
            # Integral exata de Omori:
            g_int <- 1 - (cc / (dt + cc))^(pp - 1)
            tau_clust <- sum(kappa_vals * g_int)
        } else {
            tau_clust <- 0
        }
        tau_vetor[k] <- tau_bg + tau_clust
    }
    return(tau_vetor)
}

#====Predição====

#----Previsão ETAS (curto prazo): simulação condicionada à história----
#
# Diferença em relação a simular_catalogo_etas(): ali o catálogo inteiro é
# sorteado do zero (geração 0 = eventos de fundo sorteados de mu, depois
# réplicas). Aqui, os eventos até T0 são tratados como HISTÓRIA CONHECIDA
# (fixos, não sorteados) e a simulação só gera o que acontece depois de T0:
#   (i)  réplicas futuras de cada evento histórico (mesmo que ele já tenha
#        gerado réplicas antes de T0 -- um processo de Poisson tem
#        incrementos independentes, então o número de filhos em (T0, T0+H)
#        é Poisson com média igual à integral da taxa de Omori-Utsu nesse
#        intervalo, INDEPENDENTE de quantos filhos já nasceram antes de T0);
#  (ii)  novos eventos de fundo dentro da janela [T0, T0+H];
# (iii)  toda a cascata de réplicas subsequentes desses eventos (i)+(ii),
#        usando a mesma mecânica geração-a-geração de simular_catalogo_etas.

# ---- CDF e CDF-inversa de Omori-Utsu ----
# Mesma parametrização de g_func(): g(u) = (p-1)*c^(p-1)*(u+c)^(-p), u>0.
# G(u) = integral de g de 0 a u = 1 - (c/(u+c))^(p-1).
G_omori <- function(u, cc, pp) {
    u <- pmax(u, 0)
    1 - (cc / (u + cc))^(pp - 1)
}

G_omori_inv <- function(y, cc, pp) {
    cc * ((1 - y)^(-1 / (pp - 1)) - 1)
}

# Converte um objeto catalog (pacote ETAS) ou data.frame de entrada para o
# formato padrão usado internamente: data.frame(id, time, x, y, mag):
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

# ---- Simulação de uma única realização futura, condicionada à história ----
simular_previsao_etas <- function(
                                  catalogo_historico,             # data.frame ou objeto catalog com a história observada
                                  theta,                          # Vetor nomeado: A, alpha, c, p, D, q, gamma
                                  M0      = 4.5,
                                  b       = 1.0,
                                  Mmax    = NULL,
                                  T0,                             # Tempo de origem da previsão (mesma escala do tempo do catálogo)
                                  horizonte,                      # Duração da janela de previsão [T0, T0+horizonte]
                                  mu_type = c("constant", "raster"),
                                  mu_const  = NULL,
                                  mu_raster = NULL,
                                  x_grid = NULL, y_grid = NULL,
                                  xlim = NULL, ylim = NULL,
                                  max_geracoes = 75,
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
    stopifnot(pp > 1, q > 1)
    hist <- extrair_eventos_catalog(catalogo_historico, M0)
    if (any(hist$time > T0)) {
        warning("catalogo_historico contém eventos com time > T0; Esses eventos foram descartados (não são 'história conhecida' antes da previsão).")
        hist <- hist[hist$time <= T0, , drop = FALSE]
    }
    T_fim <- T0 + horizonte
    # ---- Amostradores auxiliares ----
    amostrar_magnitudes <- function(n) {
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
    amostrar_offset_espacial <- function(mag) {
        n <- length(mag)
        if (n == 0) return(data.frame(dx = numeric(0), dy = numeric(0)))
        S_i <- D * exp(gamma * (mag - M0))
        u   <- runif(n)
        r2  <- S_i * ((1 - u)^(-1 / (q - 1)) - 1)
        r   <- sqrt(r2)
        ang <- runif(n, 0, 2 * pi)
        data.frame(dx = r * cos(ang), dy = r * sin(ang))
    }
    amostrar_local_raster <- function(n) {
        if (n == 0) return(data.frame(x = numeric(0), y = numeric(0)))
        n_x <- length(x_grid)
        probs <- as.vector(mu_raster)
        probs <- probs / sum(probs)
        idx <- sample.int(length(probs), size = n, replace = TRUE, prob = probs)
        col_idx <- ((idx - 1) %/% n_x) + 1
        row_idx <- ((idx - 1) %% n_x) + 1
        dx <- diff(x_grid)[1]; dy <- diff(y_grid)[1]
        data.frame(
            x = x_grid[row_idx] + runif(n, -dx / 2, dx / 2),
            y = y_grid[col_idx] + runif(n, -dy / 2, dy / 2)
        )
    }
    # ---- Passo 1: Réplicas futuras de eventos históricos ----
    # Para cada evento histórico i, o número de filhos que nascem em
    # (T0, T_fim) é Poisson com média kappa_i * [G(T_fim - t_i) - G(T0 - t_i)]
    # -- a diferença de CDFs garante que já não contou-se os filhos que
    # nasceram antes de T0 (esses já estão, ou deveriam estar, no próprio
    # catalogo_historico).
    n_hist <- nrow(hist)
    semente_hist <- data.frame(id = integer(0), parent_hist = character(0), geracao = integer(0),
                               time = numeric(0), x = numeric(0), y = numeric(0), mag = numeric(0),
                               origem = character(0))
    if (n_hist > 0) {
        kappa_vals <- kappa_func(hist$mag, M0, A, alpha)
        u_inf <- G_omori(T0 - hist$time, cc, pp)
        u_sup <- G_omori(T_fim - hist$time, cc, pp)
        lambda_i <- pmax(kappa_vals * (u_sup - u_inf), 0)
        n_filhos_hist <- rpois(n_hist, lambda_i)
        if (sum(n_filhos_hist) > 0) {
            idx_pais <- rep(seq_len(n_hist), times = n_filhos_hist)
            pais <- hist[idx_pais, , drop = FALSE]
            y_inf <- G_omori(T0 - pais$time, cc, pp)
            y_sup <- G_omori(T_fim - pais$time, cc, pp)
            v <- runif(nrow(pais))
            y <- y_inf + v * (y_sup - y_inf)
            u <- G_omori_inv(y, cc, pp)
            mags_filhos <- amostrar_magnitudes(nrow(pais))
            off_esp <- amostrar_offset_espacial(pais$mag)
            semente_hist <- data.frame(
                id = seq_len(nrow(pais)), parent_hist = as.character(pais$id), geracao = 1L,
                time = pais$time + u, x = pais$x + off_esp$dx, y = pais$y + off_esp$dy,
                mag = mags_filhos, origem = "descendente_historico"
            )
        }
    }
    # ---- Passo 2: Novos eventos de fundo na janela [T0, T_fim] ----
    if (mu_type == "constant") {
        stopifnot(!is.null(mu_const), !is.null(xlim), !is.null(ylim))
        area <- diff(xlim) * diff(ylim)
        n_bg <- rpois(1, mu_const * area * horizonte)
        x_bg <- runif(n_bg, xlim[1], xlim[2])
        y_bg <- runif(n_bg, ylim[1], ylim[2])
    } else {
        stopifnot(!is.null(mu_raster), !is.null(x_grid), !is.null(y_grid))
        dxg <- diff(x_grid)[1]; dyg <- diff(y_grid)[1]
        mu_total <- sum(mu_raster) * dxg * dyg
        n_bg <- rpois(1, mu_total * horizonte)
        loc <- amostrar_local_raster(n_bg)
        x_bg <- loc$x; y_bg <- loc$y
    }
    t_bg <- runif(n_bg, T0, T_fim)
    m_bg <- amostrar_magnitudes(n_bg)
    bg_futuro <- data.frame(
        id = seq_len(n_bg), parent_hist = NA_character_, geracao = 0L,
        time = t_bg, x = x_bg, y = y_bg, mag = m_bg, origem = "fundo_futuro"
    )
    # ---- Passo 3. Junta as "sementes" (geração 0 da simulação futura) ----
    sementes <- rbind(bg_futuro, semente_hist)
    if (nrow(sementes) > 0) sementes$id <- seq_len(nrow(sementes))
    proximo_id <- nrow(sementes) + 1L
    # ---- Passo 4. Propaga gerações a partir das sementes até T_fim ----
    # Mesma mecânica geração-a-geração de simular_catalogo_etas(), com T_fim
    # no lugar de T_max (os "pais" aqui já nascem dentro da janela futura,
    # então não há truncamento a considerar na descendência deles).
    catalogo_futuro <- sementes
    fila <- sementes
    g_final <- 0L
    if (nrow(fila) > 0) {
        for (g in seq_len(max_geracoes)) {
            g_final <- g
            if (nrow(fila) == 0) break
            kappa_vals <- kappa_func(fila$mag, M0, A, alpha)
            n_filhos <- rpois(nrow(fila), kappa_vals)
            if (sum(n_filhos) == 0) break
            idx_pais <- rep(seq_len(nrow(fila)), times = n_filhos)
            pais <- fila[idx_pais, , drop = FALSE]
            mags_filhos <- amostrar_magnitudes(nrow(pais))
            dt <- G_omori_inv(runif(nrow(pais)), cc, pp)
            off_esp <- amostrar_offset_espacial(pais$mag)
            filhos <- data.frame(
                id = proximo_id + seq_len(nrow(pais)) - 1L,
                parent_hist = NA_character_,
                geracao = pais$geracao + 1L,
                time = pais$time + dt,
                x = pais$x + off_esp$dx,
                y = pais$y + off_esp$dy,
                mag = mags_filhos,
                origem = "descendente_futuro"
            )
            proximo_id <- proximo_id + nrow(filhos)
            filhos <- filhos[filhos$time <= T_fim, , drop = FALSE]
            catalogo_futuro <- rbind(catalogo_futuro, filhos)
            fila <- filhos
        }
        if (g_final == max_geracoes && nrow(fila) > 0) {
            warning("Número máximo de gerações atingido na simulação de previsão -- verifique a razão de ramificação antes de confiar no resultado.")
        }
    }
    catalogo_futuro <- catalogo_futuro[order(catalogo_futuro$time), ]
    rownames(catalogo_futuro) <- NULL
    attr(catalogo_futuro, "T0") <- T0
    attr(catalogo_futuro, "T_fim") <- T_fim
    attr(catalogo_futuro, "theta") <- theta
    catalogo_futuro
}

# ---- Gera N réplicas da previsão futura (em paralelo) ----
# theta pode ser:
#  - um único vetor nomeado (previsão condicionada apenas a theta_hat, ignora
#    a incerteza de estimação -- só a estocasticidade do processo é propagada);
#  - uma lista de vetores (ex.: as ~100 estimativas do bootstrap paramétrico
#    já feito): nesse caso, a cada réplica um theta é sorteado da lista, o
#    que propaga também a incerteza de estimação de theta na previsão,
#    reaproveitando diretamente o resultado do estudo de bootstrap.
gerar_previsoes_etas <- function(
                                 catalogo_historico, theta, M0 = 4.5, b = 1.0, Mmax = NULL,
                                 T0, horizonte,
                                 mu_type = c("constant", "raster"), mu_const = NULL, mu_raster = NULL,
                                 x_grid = NULL, y_grid = NULL, xlim = NULL, ylim = NULL,
                                 n_sim = 1000, max_geracoes = 75, num_cores = 1, semente_base = NULL
                                 ) {
    mu_type <- match.arg(mu_type)
    # theta é um "ensemble" (previsão deve amostrar theta a cada réplica) se
    # for uma list() de vetores nomeados (ex.: thetas_ajustados, vindo
    # do bootstrap). Um vetor nomeado único como c(A=..., alpha=...) é
    # atomico, nao list(), entao cai no ramo "theta unico" abaixo:
    theta_eh_ensemble <- is.list(theta)
    resultados <- mclapply(seq_len(n_sim), function(s) {
        theta_s <- if (theta_eh_ensemble) theta[[sample.int(length(theta), 1)]] else theta
        semente_s <- if (!is.null(semente_base)) semente_base + s else NULL
        simular_previsao_etas(
            catalogo_historico = catalogo_historico,
            theta = theta_s, M0 = M0, b = b, Mmax = Mmax,
            T0 = T0, horizonte = horizonte,
            mu_type = mu_type, mu_const = mu_const, mu_raster = mu_raster,
            x_grid = x_grid, y_grid = y_grid, xlim = xlim, ylim = ylim,
            max_geracoes = max_geracoes, checar_criticidade = FALSE,
            semente = semente_s
        )
    }, mc.cores = num_cores)
    list(simulacoes = resultados, T0 = T0, horizonte = horizonte, n_sim = n_sim)
}

# ---- Agregação temporal: contagem esperada de eventos por bin de tempo ----
agregar_previsao_temporal <- function(previsoes, largura_bin = 1) {
    T0 <- previsoes$T0
    horizonte <- previsoes$horizonte
    lista_futuros <- previsoes$simulacoes
    quebras <- seq(T0, T0 + horizonte, by = largura_bin)
    if (max(quebras) < T0 + horizonte) quebras <- c(quebras, T0 + horizonte)
    n_bins <- length(quebras) - 1
    contagens <- vapply(lista_futuros, function(cat_fut) {
        if (nrow(cat_fut) == 0) return(rep(0, n_bins))
        as.numeric(table(cut(cat_fut$time, breaks = quebras, include.lowest = TRUE, right = FALSE)))
    }, numeric(n_bins))
    # 'contagens' é uma matriz n_bins x n_sim
    data.frame(
        bin_inicio = quebras[-length(quebras)],
        bin_fim    = quebras[-1],
        media      = rowMeans(contagens),
        mediana    = apply(contagens, 1, median),
        p025       = apply(contagens, 1, quantile, probs = 0.025),
        p975       = apply(contagens, 1, quantile, probs = 0.975)
    )
}

# ---- Agregação espacial: taxa esperada e P(n_eventos >=1 evento) por célula do grid ----
agregar_previsao_espacial <- function(previsoes, x_grid, y_grid) {
    lista_futuros <- previsoes$simulacoes
    n_sim <- length(lista_futuros)
    n_x <- length(x_grid); n_y <- length(y_grid)
    dx <- diff(x_grid)[1]; dy <- diff(y_grid)[1]
    quebras_x <- c(x_grid - dx / 2, max(x_grid) + dx / 2)
    quebras_y <- c(y_grid - dy / 2, max(y_grid) + dy / 2)
    soma_contagens <- matrix(0, nrow = n_x, ncol = n_y)
    for (cat_fut in lista_futuros) {
        if (nrow(cat_fut) == 0) next
        ix <- cut(cat_fut$x, breaks = quebras_x, include.lowest = TRUE, labels = FALSE)
        iy <- cut(cat_fut$y, breaks = quebras_y, include.lowest = TRUE, labels = FALSE)
        validos <- !is.na(ix) & !is.na(iy)
        if (any(validos)) {
            tab <- table(factor(ix[validos], levels = seq_len(n_x)), factor(iy[validos], levels = seq_len(n_y)))
            soma_contagens <- soma_contagens + as.matrix(tab)
        }
    }
    contagem_media <- soma_contagens / n_sim
    list(
        x_grid = x_grid,
        y_grid = y_grid,
        contagem_media    = contagem_media,   # Número esperado de eventos por célula
        prob_pelo_menos_1 = 1 - exp(-contagem_media)  # Aproximacao de Poisson
    )
}

#====Perfil de verossimilhança para um modelo já ajustado (ajustar_etas_proprio)====

#----Versão RÁPIDA: fatia da verossimilhança (E-step fixo no ajuste final)----
# Custo: 1 chamada de calcular_probabilidades() + |grade_valores| avaliações
# de log_lik_etas() (baratas). Não reotimiza os outros parâmetros -- então
# tende a EXAGERAR o quão rápido a verossimilhança piora ao afastar do MLE,
# já que ignora a compensação de outros parâmetros correlacionados (ex.: D e
# q compensando gamma, como vimos no estudo de sensibilidade). Útil como
# primeira checagem visual rápida, não para inferência formal.
perfilar_rapido <- function(ajuste, df_eventos, x_grid, y_grid, M0 = 4.5, num_cores = 1,
                            parametro, grade_valores) {
    passo_e <- calcular_probabilidades(df_eventos, ajuste$theta, M0, ajuste$raster, x_grid, y_grid, num_cores)
    logvero <- vapply(grade_valores, function(v) {
        theta_v <- ajuste$theta
        theta_v[parametro] <- v
        -log_lik_etas(theta_v, df_eventos, passo_e$probabilidades, passo_e$prob_total,
                      ajuste$raster, x_grid, y_grid, M0 = M0)
    }, numeric(1))
    data.frame(valor = grade_valores, log_verossimilhanca = logvero)
}

#---- Versão COMPLETA: perfil de verossimilhança de verdade----
# Para cada valor de 'parametro' na grade: fixa esse parâmetro, recalcula o
# passo E consistente com ele (1 chamada de calcular_probabilidades por
# ponto da grade -- não a cada iteração interna do optim, para manter o
# custo controlado), e reotimiza os 6 parâmetros restantes via optim() com
# gradiente analítico. O raster de fundo é mantido fixo no valor do ajuste
# original (reestimá-lo também a cada ponto da grade seria o "EM completo
# refeito por ponto" -- mais fiel, mas várias vezes mais caro; ver nota
# abaixo). Os pontos da grade são independentes entre si -- paralelizável
# via mc.cores.
perfilar_completo <- function(ajuste, df_eventos, x_grid, y_grid, M0 = 4.5, num_cores_e = 1,
                              parametro, grade_valores,
                              lower_full = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-6),
                              upper_full = c(15, 10, 15, 15, 15, 15, 10),
                              mc_cores_perfil = 1) {
    nomes_theta <- names(ajuste$theta)
    idx_fixo    <- which(nomes_theta == parametro)
    idx_livre   <- setdiff(seq_along(nomes_theta), idx_fixo)
    r_fixo      <- ajuste$raster
    names(lower_full) <- names(upper_full) <- nomes_theta
    resultado <- mclapply(grade_valores, function(v) {
        theta_ponto <- ajuste$theta
        theta_ponto[idx_fixo] <- v
        # Passo E consistente com esse valor de 'parametro' (uma vez por
        # ponto da grade, não a cada avaliação interna do otimizador):
        passo_e <- calcular_probabilidades(df_eventos, theta_ponto, M0, r_fixo, x_grid, y_grid, num_cores_e)
        obj <- function(theta_livre) {
            theta_completo <- theta_ponto
            theta_completo[idx_livre] <- theta_livre
            log_lik_etas(theta_completo, df_eventos, passo_e$probabilidades, passo_e$prob_total,
                        r_fixo, x_grid, y_grid, M0 = M0)
        }
        grad_obj <- function(theta_livre) {
            theta_completo <- theta_ponto
            theta_completo[idx_livre] <- theta_livre
            gradiente_etas_completo(theta_completo, df_eventos, passo_e$probabilidades, passo_e$prob_total, M0 = M0)[idx_livre]
        }
        otim <- tryCatch(
            optim(theta_ponto[idx_livre], obj, gr = grad_obj, method = "L-BFGS-B",
                 lower = lower_full[idx_livre], upper = upper_full[idx_livre],
                 control = list(maxit = 500)),
            error = function(e) list(value = NA, convergence = NA_integer_)
        )
        list(valor = v, log_verossimilhanca = if (!is.na(otim$value)) -otim$value else NA,
            convergence = otim$convergence)
    }, mc.cores = mc_cores_perfil)
    data.frame(
        valor               = vapply(resultado, function(x) x$valor, numeric(1)),
        log_verossimilhanca = vapply(resultado, function(x) x$log_verossimilhanca, numeric(1)),
        convergence         = vapply(resultado, function(x) x$convergence, integer(1))
    )
}
