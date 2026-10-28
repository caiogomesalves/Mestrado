library(dplyr)
library(lubridate)
library(parallel)
library(tools)
library(ETAS)

#====ETAS::etas() com salvaguarda de tempo-limite (processo isolado)====
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

# ETAS::etas() roda em código C/C++ compilado (portado de Fortran) sem trava
# de caixa no espaço de parâmetros. Se divergir, pode ficar preso sem nunca
# checar interrupção do R -- por isso a única salvaguarda confiável é isolar
# a chamada em um processo filho e matá-lo a nível de SO se não retornar a
# tempo.
rodar_etas_seguro <- function(catalogo_etas, param0, tempo_limite_seg = 300, ...) {
    job <- mcparallel({
        tryCatch(
            etas(catalogo_etas, param0 = param0, verbose = FALSE, plot.it = FALSE, nthreads = 15, ...),
            error = function(e) list(erro_r = conditionMessage(e))
        )
    }, silent = TRUE)
    resultado <- mccollect(job, wait = TRUE, timeout = tempo_limite_seg)
    if (is.null(resultado)) {
        pskill(job$pid, signal = SIGKILL)
        mccollect(job)
        return(list(status = "timeout", ajuste = NULL,
                    mensagem = sprintf("Processo morto após %ds sem retornar (parâmetros provavelmente divergiram).", tempo_limite_seg)))
    }
    ajuste <- resultado[[1]]
    if (is.list(ajuste) && !is.null(ajuste$erro_r)) {
        return(list(status = "erro_r", ajuste = NULL, mensagem = ajuste$erro_r))
    }
    if (inherits(ajuste, "try-error")) {
        return(list(status = "erro_r", ajuste = NULL, mensagem = paste(ajuste)))
    }
    params_finais <- ajuste$param
    if (any(!is.finite(params_finais)) || any(abs(params_finais) > 1e6)) {
        return(list(status = "divergiu", ajuste = ajuste,
                    mensagem = "Parâmetros finais não-finitos ou implausivelmente grandes."))
    }
    list(status = "ok", ajuste = ajuste, mensagem = NULL)
}

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

gerar_lhs <- function(n, k, semente = NULL) {
    if (!is.null(semente)) set.seed(semente)
    U <- matrix(runif(n * k), nrow = n, ncol = k)
    for (j in seq_len(k)) {
        perm <- sample.int(n)
        U[, j] <- (perm - 1 + U[, j]) / n
    }
    U
}

#====Loop de comparação: método próprio vs. ETAS::etas(), mesmos pontos de partida====

theta_verdadeiro <- c("A" = 0.5, "alpha" = 1.5, "c" = 0.25,
                     "p" = 1.5, "D" = 0.1, "q" = 2, "gamma" = 0.1)

nomes_theta <- names(theta_verdadeiro)

# ---- Design 1: "um fator de cada vez" (OAT) -- igual ao original, mas
# corrigido para respeitar p,q > 1 ----
mult_grid_oat <- seq(0.25, 2, length.out = 10)

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

n_multifator <- 50
faixa_mult    <- c(0.25, 4)

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

#----Ajuste usando ETAS::etas----

resultados_etas_pkg <- vector("list", n_pontos_partida)

mu_inicial               <- 0.1   # chute inicial p/ taxa de fundo constante do ETAS::etas() (ajustar)
tempo_limite_por_ajuste  <- 600   # segundos -- calibrar rodando 1-2 ajustes "bons" primeiro

for (k in 1:n_pontos_partida) {
    theta0_k <- theta_matrix_design[k, ]
    # ---- ETAS::etas(), com salvaguarda de tempo-limite ----
    # ETAS::etas() espera um vetor com TODOS os parâmetros, incluindo mu
    # (taxa de fundo constante) -- confira a ordem exata em ?ETAS::etas
    # para a versão instalada antes de rodar; ajuste se for diferente.
    param0_completo <- c(mu = mu_inicial, theta0_k)
    resultados_etas_pkg[[k]] <- rodar_etas_seguro(
        cat_temp_iran, param0 = param0_completo,
        tempo_limite_seg = tempo_limite_por_ajuste, no.itr = 11
    )
    cat("==========Iteração", k, "completa==========")
    cat(sprintf("\n"))
    if (k %% 10 == 0 || k == n_pontos_partida) {
        save(resultados_etas_pkg, k,
             file = "checkpoint_comparacao_etas.RData", compress = TRUE)
    }
}

status_etas_pkg <- vapply(resultados_etas_pkg, function(x) x$status, character(1))
cat("\n---- Resumo dos status do ETAS::etas() por ponto de partida ----\n")
print(table(status_etas_pkg))
