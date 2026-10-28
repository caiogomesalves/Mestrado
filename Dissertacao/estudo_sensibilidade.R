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

region_iran <- list(lat = c(26, 25, 29, 38, 35),
                    long = c(52, 59, 58, 45, 43))

# Catálogo:
iran.cat <- catalog(iran.quakes,
                    study.start = "1991/01/01", study.end = "2011/01/01",
                    region.poly = region_iran, mag.threshold = 4.5)

# Chute inicial parâmetros:
theta_iran <- c("A" = 0.23, "alpha" = 2.8, "c" = 0.022,
                "p" = 1.12, "D" = 0.012, "q" = 2.3, "gamma" = 0.03)

thetas_iran <- list()
thetas_iran[[1]] <- theta_iran

# Número total de eventos:
n_events_iran <- nrow(iran.cat$revents)

# Definir as sequências do grid explicitamente
x_grid_iran <- seq(min(iran.cat$longlat.coord$long), max(iran.cat$longlat.coord$long), length.out = 128)
y_grid_iran <- seq(min(iran.cat$longlat.coord$lat), max(iran.cat$longlat.coord$lat), length.out = 128)

# Raster inicial:
r_iran <- matrix(1, nrow = 128, ncol = 128)
rel_clust_iran <- matrix(0, nrow = 128, ncol = 128)

# Normaliza o raster para a integral ser 1:
# r_iran <- r_iran/(sum(r_iran))

# Na lista de rasters, guardaremos as matrizes:
rasters_iran <- list()
rast_rel_clust_iran <- list()
rasters_iran[[1]] <- r_iran

hessianas_iran <- list()

# Data.frame com os eventos para funções:
df_temp_iran <- data.frame(x = iran.cat$longlat.coord$long,
                           y = iran.cat$longlat.coord$lat,
                           mag = iran.cat$revents[, "mm"] + 4.5,
                           time = iran.cat$revents[, "tt"])

# Distâncias para kernel adaptativo:
d_j_iran <- d_j(iran.cat$longlat.coord[, c("long", "lat")])

num_cores <- 15

# Cálculo (paralelo) das probabilidades p_ij e das probabilidades totais:
passo_e_iran <- calcular_probabilidades(df_temp_iran, theta_iran, 4.5, r_iran, x_grid_iran, y_grid_iran, num_cores)
probabilidades_iran <- passo_e_iran$probabilidades
prob_total_iran <- passo_e_iran$prob_total

# Estima o background rate:
r_iran <- estimate_background_kernel(x_grid_iran,
                                     y_grid_iran,
                                     data.frame(x = iran.cat$longlat.coord$long,
                                                y = iran.cat$longlat.coord$lat,
                                                rho = prob_total_iran),
                                     d_j_iran,
                                     T_total = max(df_temp_iran$time))

rel_clust_iran <- estimate_relative_clustering(x_grid_iran,
                                               y_grid_iran,
                                               data.frame(x = iran.cat$longlat.coord$long,
                                                          y = iran.cat$longlat.coord$lat,
                                                          rho = prob_total_iran),
                                               d_j_iran)$C

# Otimiza os parâmetros theta:
otimizacao_iran <- optim(theta_iran,
                         log_lik_etas,
                         gr = gradiente_etas_completo,
                         df_eventos = df_temp_iran,
                         matriz_p_ij = probabilidades_iran,
                         prob_total = prob_total_iran,
                         raster_mu = r_iran, x_grid = x_grid_iran, y_grid = y_grid_iran,
                         method = "L-BFGS-B",
                         lower = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-10),
                         upper = c(15, 10, 15, 15, 15, 15, 10),
                         control = list(maxit = 1000),
                         hessian = T)

theta_iran <- otimizacao_iran$par
thetas_iran[[2]] <- theta_iran

hessianas_iran[[1]] <- otimizacao_iran$hessian

rast_rel_clust_iran[[1]] <- rel_clust_iran

rasters_iran[[2]] <- r_iran

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
while (!convergiu_em(rasters_iran[[i]], rasters_iran[[i - 1]], tol_rel_raster, tol_abs_raster) |
       !convergiu_em(thetas_iran[[i]], thetas_iran[[i - 1]], tol_rel_theta, tol_abs_theta)){
    if (i >= max_iter_em) {
        warning(sprintf("EM (ajuste real) atingiu max_iter_em = %d sem convergir dentro da tolerancia -- interrompendo.", max_iter_em))
        break
    }
    # Calcula pares de probabilidades (passo E):
    passo_e_iran <- calcular_probabilidades(df_temp_iran, theta_iran, 4.5, r_iran, x_grid_iran, y_grid_iran, num_cores)
    probabilidades_iran <- passo_e_iran$probabilidades
    prob_total_iran <- passo_e_iran$prob_total
    # Estima parâmetros theta:
    otimizacao_iran <- optim(theta_iran,
                             log_lik_etas,
                             gr = gradiente_etas_completo,
                             df_eventos = df_temp_iran,
                             matriz_p_ij = probabilidades_iran,
                             prob_total = prob_total_iran,
                             raster_mu = r_iran, x_grid = x_grid_iran, y_grid = y_grid_iran,
                             method = "L-BFGS-B",
                             lower = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-10),
                             upper = c(15, 10, 15, 15, 15, 15, 10),
                             control = list(maxit = 1000),
                             hessian = T)
    if (otimizacao_iran$convergence != 0) {
        warning(sprintf("optim() nao convergiu na iteracao %d do ajuste real (codigo = %d): %s",
                        i, otimizacao_iran$convergence, otimizacao_iran$message))
    }
    theta_iran <- otimizacao_iran$par
    thetas_iran[[i + 1]] <- theta_iran
    hessianas_iran[[i]] <- otimizacao_iran$hessian
    # Estima background rate:
    r_iran <- estimate_background_kernel(x_grid_iran,
                                         y_grid_iran,
                                         data.frame(x = iran.cat$longlat.coord$long,
                                                    y = iran.cat$longlat.coord$lat,
                                                    rho = prob_total_iran),
                                         d_j_iran,
                                         T_total = max(df_temp_iran$time))
    rasters_iran[[i + 1]] <- r_iran
    # Estima clustering relativo:
    rel_clust_iran <- estimate_relative_clustering(x_grid_iran,
                                                   y_grid_iran,
                                                   data.frame(x = iran.cat$longlat.coord$long,
                                                              y = iran.cat$longlat.coord$lat,
                                                              rho = prob_total_iran),
                                                   d_j_iran)$C
    rast_rel_clust_iran[[i]] <- rel_clust_iran
    # Mostra progresso do loop:
    print(i + 1)
    # Atualiza índice:
    i <- i + 1
}

#====Estudo de sensibilidade ao chute inicial (Irã)====
# Simula UM catálogo fixo a partir de um theta "verdadeiro" conhecido, e
# reajusta esse mesmo catálogo a partir de vários pontos de partida
# diferentes para o EM/optim(). A pergunta respondida é: "partindo de pontos
# de partida diferentes, o ajuste converge de volta para a verdade, ou fica
# preso perto do chute inicial?" -- isto é DIFERENTE de um bootstrap (que usa
# muitos catálogos simulados independentes para medir erro amostral); aqui
# há só um catálogo, e a variabilidade entre os ajustes finais reflete
# sensibilidade NUMÉRICA do otimizador ao ponto de partida, não incerteza
# estatística. Por isso os nomes de variáveis abaixo evitam o vocabulário
# de "bootstrap"/"réplica" usado no estudo de erro padrão.

# ---- Theta "verdadeiro" usado para simular os dados (definido uma única
# vez e reutilizado em toda a análise -- corrige o bug em que a checagem de
# criticidade usava theta_final_iran, do ajuste no catálogo REAL, em vez do
# theta que de fato gerou o catálogo sintético testado aqui) ----
theta_verdadeiro <- c("A" = 0.5, "alpha" = 1.5, "c" = 0.25,
                     "p" = 1.5, "D" = 0.1, "q" = 2, "gamma" = 0.1)

nomes_theta <- names(theta_verdadeiro)

raster_final_iran <- rasters_iran[[length(rasters_iran)]]

# Checagem de criticidade do modelo que efetivamente gera os dados simulados
# (agora aponta para theta_verdadeiro, não mais para theta_final_iran):
verificar_criticidade(theta_verdadeiro, b = 2.0)

# Grid espacial (fixo, calculado uma única vez):
x_grid_iran <- seq(min(iran.cat$longlat.coord$long), max(iran.cat$longlat.coord$long), length.out = 128)
y_grid_iran <- seq(min(iran.cat$longlat.coord$lat), max(iran.cat$longlat.coord$lat), length.out = 128)

# ---- Perturbação de theta respeitando as restrições do modelo ----
# Para A, alpha, c, D, gamma (que só precisam ser > 0): valor = verdadeiro * multiplicador.
# Para p e q (que precisam ser > 1 -- restrição rígida da lei de Omori-Utsu e
# do kernel espacial de Pareto): perturbamos o EXCESSO acima de 1, não o
# valor bruto, para nunca cruzar essa fronteira:
#     valor = 1 + (verdadeiro - 1) * multiplicador
# (o design OAT original multiplicava p e q diretamente por 0.25, gerando
# p0 = 0.375 e q0 = 0.5 -- ambos abaixo do mínimo permitido de 1. Esse era
# um bug real do script enviado, corrigido aqui.)
perturbar_theta <- function(theta_verdadeiro, multiplicadores) {
    stopifnot(all(names(multiplicadores) == names(theta_verdadeiro)))
    theta_novo <- theta_verdadeiro
    params_excesso <- c("p", "q")
    for (nome in names(theta_verdadeiro)) {
        if (nome %in% params_excesso) {
            theta_novo[nome] <- 1 + (theta_verdadeiro[nome] - 1) * multiplicadores[nome]
        } else {
            theta_novo[nome] <- theta_verdadeiro[nome] * multiplicadores[nome]
        }
    }
    theta_novo
}

# ---- Design 1: "um fator de cada vez" (OAT) -- igual ao original, mas
# corrigido para respeitar p,q > 1 ----
mult_grid_oat <- seq(0.25, 2, length.out = 10)  # multiplicadores testados por parâmetro

theta_matrix_oat <- matrix(nrow = 0, ncol = length(nomes_theta), dimnames = list(NULL, nomes_theta))
metadata_oat <- data.frame(tipo = character(0), parametro_variado = character(0), multiplicador = numeric(0))

for (par_variado in nomes_theta) {
    for (m in mult_grid_oat) {
        mult <- setNames(rep(1, length(nomes_theta)), nomes_theta)
        mult[par_variado] <- m
        theta0 <- perturbar_theta(theta_verdadeiro, mult)
        theta_matrix_oat <- rbind(theta_matrix_oat, theta0)
        metadata_oat <- rbind(metadata_oat, data.frame(tipo = "OAT", parametro_variado = par_variado, multiplicador = m))
    }
}
rownames(theta_matrix_oat) <- NULL

# ---- Design 2: perturbação de MÚLTIPLOS fatores simultaneamente ----
# Em vez de fixar 6 parâmetros na verdade e variar só 1 (OAT), aqui sorteamos
# um multiplicador para CADA um dos 7 parâmetros ao mesmo tempo, usando
# Latin Hypercube Sampling (LHS) em escala log: para n pontos de partida e k
# parâmetros, cada uma das k dimensões é dividida em n faixas, e cada faixa é
# usada exatamente uma vez (uma permutação aleatória por dimensão) -- isso
# garante cobertura equilibrada do espaço de 7 dimensões com muito menos
# pontos do que um fatorial completo (2^7 = 128 combinações só para testar
# "baixo/alto" em cada parâmetro, sem nem cobrir o interior do espaço).
# O objetivo é testar interações: será que o ajuste continua robusto quando
# VÁRIOS parâmetros começam errados ao mesmo tempo (ex.: A e alpha, que são
# historicamente correlacionados na verossimilhança do ETAS), não só quando
# um erra por vez.
gerar_lhs <- function(n, k, semente = NULL) {
    if (!is.null(semente)) set.seed(semente)
    U <- matrix(runif(n * k), nrow = n, ncol = k)
    for (j in seq_len(k)) {
        perm <- sample.int(n)
        U[, j] <- (perm - 1 + U[, j]) / n
    }
    U
}

n_multifator <- 50
faixa_mult    <- c(0.25, 4)  # multiplicador minimo/maximo testado, em escala log

lhs_unit <- gerar_lhs(n_multifator, length(nomes_theta), semente = 321)
colnames(lhs_unit) <- nomes_theta

log_min <- log(faixa_mult[1]); log_max <- log(faixa_mult[2])
mult_lhs <- exp(log_min + lhs_unit * (log_max - log_min))  # multiplicadores em escala log (equilibra sub/superestimar)

theta_matrix_lhs <- t(apply(mult_lhs, 1, function(mult_linha) {
    perturbar_theta(theta_verdadeiro, setNames(mult_linha, nomes_theta))
}))
colnames(theta_matrix_lhs) <- nomes_theta
rownames(theta_matrix_lhs) <- NULL

metadata_lhs <- data.frame(tipo = rep("LHS_multifator", n_multifator),
                           parametro_variado = rep(NA_character_, n_multifator),
                           multiplicador = rep(NA_real_, n_multifator))

# ---- Combina os dois designs em uma única matriz de pontos de partida ----
theta_matrix_design <- rbind(theta_matrix_oat, theta_matrix_lhs)
metadata_design      <- rbind(metadata_oat, metadata_lhs)
n_pontos_partida     <- nrow(theta_matrix_design)

# ---- Simula UM catálogo fixo a partir de theta_verdadeiro ----
cat_sim_iran <- simular_catalogo_etas(
    theta   = theta_verdadeiro,
    M0      = 4.5,
    b       = 1.5,
    Mmax    = 10,
    T_max   = 3650,
    mu_type = "raster",
    mu_raster = raster_final_iran,
    x_grid = x_grid_iran,
    y_grid = y_grid_iran,
    checar_criticidade = FALSE,
    semente = 123
)

cat_temp_iran <- catalog(
    cat_sim_iran %>%
    mutate(t = ddays(time),
           data = POSIXct(1) + t,
           long = x, lat = y,
           date = date(data),
           time = format(data, format = "%H:%M:%S")),
    region.poly = iran.cat$region.poly
)

df_temp_iran_temp <- data.frame(x = cat_temp_iran$longlat.coord$long,
                                y = cat_temp_iran$longlat.coord$lat,
                                mag = cat_temp_iran$revents[, "mm"] + 4.5,
                                time = cat_temp_iran$revents[, "tt"])
df_temp_iran_temp$time <- df_temp_iran_temp$time - min(df_temp_iran_temp$time)

# O catálogo é o MESMO em todos os pontos de partida -- diferente do
# bootstrap (onde cada réplica tinha dados novos), aqui d_j e o número de
# eventos não mudam entre pontos de partida, então são calculados uma única
# vez fora do loop (no script original, isso era recalculado a cada uma das
# 70 iterações sem necessidade):
d_j_iran_temp       <- d_j(cat_temp_iran$longlat.coord[, c("long", "lat")])
n_events_iran_simu  <- nrow(cat_temp_iran$revents)

num_cores         <- 15
max_iter_em       <- 100
checkpoint_a_cada <- 10

thetas_finais_sens      <- vector("list", n_pontos_partida)
backgrounds_finais_sens <- vector("list", n_pontos_partida)
convergencia_sens       <- vector("list", n_pontos_partida)
falhas_sens             <- vector("list", n_pontos_partida)

for (k in 1:n_pontos_partida) {
    resultado_k <- tryCatch({
        theta_iran_simu <- theta_matrix_design[k, ]
        thetas_iran_temp <- list()
        thetas_iran_temp[[1]] <- theta_iran_simu
        r_iran_temp <- matrix(1, nrow = 128, ncol = 128)
        rasters_iran_temp <- list()
        rasters_iran_temp[[1]] <- r_iran_temp

        passo_e_temp <- calcular_probabilidades(df_temp_iran_temp, theta_iran_simu, 4.5, r_iran_temp, x_grid_iran, y_grid_iran, num_cores)
        probabilidades_iran_temp <- passo_e_temp$probabilidades
        prob_total_iran_temp <- passo_e_temp$prob_total

        r_iran_temp <- estimate_background_kernel(x_grid_iran,
                                                  y_grid_iran,
                                                  data.frame(x = cat_temp_iran$longlat.coord$long,
                                                             y = cat_temp_iran$longlat.coord$lat,
                                                             rho = prob_total_iran_temp),
                                                  d_j_iran_temp,
                                                  T_total = max(df_temp_iran_temp$time))

        otimizacao_iran_temp <- optim(theta_iran_simu,
                                      log_lik_etas,
                                      gr = gradiente_etas_completo,
                                      df_eventos = df_temp_iran_temp,
                                      matriz_p_ij = probabilidades_iran_temp,
                                      prob_total = prob_total_iran_temp,
                                      raster_mu = r_iran_temp, x_grid = x_grid_iran, y_grid = y_grid_iran,
                                      method = "L-BFGS-B",
                                      lower = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-6),
                                      upper = c(15, 10, 15, 15, 15, 15, 10),
                                      control = list(maxit = 1000),
                                      hessian = T)
        conv_codes <- otimizacao_iran_temp$convergence
        theta_iran_simu <- otimizacao_iran_temp$par
        thetas_iran_temp[[2]] <- theta_iran_simu
        rasters_iran_temp[[2]] <- r_iran_temp

        i_em <- 2
        while (!convergiu_em(rasters_iran_temp[[i_em]], rasters_iran_temp[[i_em - 1]], tol_rel_raster, tol_abs_raster) |
               !convergiu_em(thetas_iran_temp[[i_em]], thetas_iran_temp[[i_em - 1]], tol_rel_theta, tol_abs_theta)) {

            if (i_em >= max_iter_em) {
                warning(sprintf("Ponto de partida %d: EM atingiu max_iter_em = %d sem convergir -- interrompendo.", k, max_iter_em))
                break
            }

            passo_e_temp <- calcular_probabilidades(df_temp_iran_temp, theta_iran_simu, 4.5, r_iran_temp, x_grid_iran, y_grid_iran, num_cores)
            probabilidades_iran_temp <- passo_e_temp$probabilidades
            prob_total_iran_temp <- passo_e_temp$prob_total

            otimizacao_iran_temp <- optim(theta_iran_simu,
                                          log_lik_etas,
                                          gr = gradiente_etas_completo,
                                          df_eventos = df_temp_iran_temp,
                                          matriz_p_ij = probabilidades_iran_temp,
                                          prob_total = prob_total_iran_temp,
                                          raster_mu = r_iran_temp, x_grid = x_grid_iran, y_grid = y_grid_iran,
                                          method = "L-BFGS-B",
                                          lower = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-6),
                                          upper = c(15, 10, 15, 15, 15, 15, 10),
                                          control = list(maxit = 1000),
                                          hessian = T)
            conv_codes <- c(conv_codes, otimizacao_iran_temp$convergence)
            theta_iran_simu <- otimizacao_iran_temp$par
            thetas_iran_temp[[i_em + 1]] <- theta_iran_simu

            r_iran_temp <- estimate_background_kernel(x_grid_iran,
                                                      y_grid_iran,
                                                      data.frame(x = cat_temp_iran$longlat.coord$long,
                                                                 y = cat_temp_iran$longlat.coord$lat,
                                                                 rho = prob_total_iran_temp),
                                                      d_j_iran_temp,
                                                      T_total = max(df_temp_iran_temp$time))
            rasters_iran_temp[[i_em + 1]] <- r_iran_temp
            cat(sprintf("[Ponto de partida %d] iteracao EM %d\n", k, i_em + 1))
            i_em <- i_em + 1
        }

        list(
            theta = thetas_iran_temp[[length(thetas_iran_temp)]],
            raster = rasters_iran_temp[[length(rasters_iran_temp)]],
            n_iter_em = i_em,
            convergence_codes = conv_codes,
            convergiu = all(conv_codes == 0) && i_em < max_iter_em
        )

    }, error = function(e) {
        message(sprintf("Ponto de partida %d falhou: %s", k, conditionMessage(e)))
        list(erro = conditionMessage(e))
    })

    if (!is.null(resultado_k$erro)) {
        falhas_sens[[k]] <- resultado_k$erro
        thetas_finais_sens[[k]] <- NA
        backgrounds_finais_sens[[k]] <- NA
        convergencia_sens[[k]] <- NA
    } else {
        thetas_finais_sens[[k]] <- resultado_k$theta
        backgrounds_finais_sens[[k]] <- resultado_k$raster
        convergencia_sens[[k]] <- list(n_iter_em = resultado_k$n_iter_em,
                                       convergence_codes = resultado_k$convergence_codes,
                                       convergiu = resultado_k$convergiu)
    }

    cat(sprintf("==== Ponto de partida %d/%d concluído ====\n", k, n_pontos_partida))

    if (k %% checkpoint_a_cada == 0 || k == n_pontos_partida) {
        save(thetas_finais_sens, backgrounds_finais_sens, convergencia_sens, falhas_sens, k,
             file = "checkpoint_sensibilidade_iran.RData", compress = TRUE)
    }
}

#----Análise final: o ajuste recupera o theta verdadeiro? ----
# Esse é o passo que faltava no script original: comparar cada ajuste final
# com theta_verdadeiro (não só armazenar os ajustes), e estratificar por qual
# parâmetro foi perturbado (design OAT) ou pelo tipo de design (multifator).

idx_validos <- which(vapply(convergencia_sens, function(x) !identical(x, NA) && isTRUE(x$convergiu), logical(1)))

if (length(idx_validos) < n_pontos_partida) {
    warning(sprintf("%d de %d pontos de partida foram descartados (erro no ajuste ou EM sem convergência).",
                    n_pontos_partida - length(idx_validos), n_pontos_partida))
}

theta_matrix_final <- do.call(rbind, thetas_finais_sens[idx_validos])
metadata_validos   <- metadata_design[idx_validos, ]

# Erro relativo de cada ajuste final em relação ao theta verdadeiro
# (positivo = superestimou, negativo = subestimou):
theta_verdadeiro_rep <- matrix(theta_verdadeiro[colnames(theta_matrix_final)],
                               nrow = nrow(theta_matrix_final), ncol = ncol(theta_matrix_final), byrow = TRUE)
erro_relativo <- (theta_matrix_final - theta_verdadeiro_rep) / theta_verdadeiro_rep
colnames(erro_relativo) <- colnames(theta_matrix_final)

resumo_sensibilidade <- data.frame(
    tipo              = metadata_validos$tipo,
    parametro_variado = metadata_validos$parametro_variado,
    multiplicador     = metadata_validos$multiplicador,
    erro_abs_medio    = rowMeans(abs(erro_relativo))   # média do |erro relativo| entre os 7 parâmetros
)

# Resumo por parâmetro perturbado (design OAT) -- responde "o ajuste é mais
# sensível ao chute inicial de qual parâmetro?":
resumo_oat <- aggregate(erro_abs_medio ~ parametro_variado,
                        data = subset(resumo_sensibilidade, tipo == "OAT"),
                        FUN = function(x) c(media = mean(x), maximo = max(x)))

# Resumo do design multifator -- responde "o ajuste continua robusto quando
# vários parâmetros começam errados ao mesmo tempo?":
resumo_lhs <- subset(resumo_sensibilidade, tipo == "LHS_multifator")

# Pontos de partida "problemáticos": erro médio acima de 20% mesmo após
# convergência -- candidatos a ótimo local, não só ruído numérico:
limiar_erro_problematico <- 0.20
pontos_problematicos <- resumo_sensibilidade[resumo_sensibilidade$erro_abs_medio > limiar_erro_problematico, ]

cat("\n---- Erro absoluto médio por parâmetro perturbado (design OAT) ----\n")
print(resumo_oat)
cat("\n---- Erro absoluto médio (design multifator/LHS) ----\n")
cat(sprintf("Média: %.4f | Mediana: %.4f | Máximo: %.4f\n",
           mean(resumo_lhs$erro_abs_medio), median(resumo_lhs$erro_abs_medio), max(resumo_lhs$erro_abs_medio)))
cat(sprintf("\n%d de %d pontos de partida válidos tiveram erro médio acima de %.0f%% -- possíveis ótimos locais.\n",
           nrow(pontos_problematicos), length(idx_validos), limiar_erro_problematico * 100))

#----Salvando os Resultados----

save(
    theta_verdadeiro,
    theta_matrix_design,
    metadata_design,
    thetas_finais_sens,
    backgrounds_finais_sens,
    convergencia_sens,
    falhas_sens,
    theta_matrix_final,
    erro_relativo,
    resumo_sensibilidade,
    resumo_oat,
    resumo_lhs,
    pontos_problematicos,
    file = "resultados_sensibilidade_iran.RData",
    compress = TRUE
)

print("Estudo de sensibilidade finalizado com sucesso. Arquivos salvos.")
