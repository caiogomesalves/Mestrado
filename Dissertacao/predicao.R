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
#
# Observação de implementação: os amostradores auxiliares (magnitude, offset
# espacial, localização de fundo via raster) estão duplicados aqui de forma
# autocontida, em vez de reaproveitar os closures internos de
# simular_catalogo_etas(). Isso é proposital: evita qualquer risco de alterar
# a ordem dos sorteios aleatórios dentro da função já usada no estudo de
# bootstrap (o que mudaria resultados reproduzidos com uma mesma semente).

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
# formato padrão usado internamente: data.frame(id, time, x, y, mag).
# - Se for um objeto "catalog", extrai os eventos com flag em {0,1} (mesma
#   convenção usada no restante do script) e reconstrói a magnitude absoluta
#   somando M0 (no pacote ETAS, "mm" é a magnitude relativa ao limiar).
# - Se for um data.frame, exige as colunas time, x, y, mag (id é opcional;
#   se ausente, é criado a partir da ordem das linhas).
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
            stop(sprintf("catalogo_historico esta sem a(s) coluna(s): %s", paste(faltando, collapse = ", ")))
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
                                  catalogo_historico,             # data.frame ou objeto catalog com a historia observada
                                  theta,                          # vetor nomeado: A, alpha, c, p, D, q, gamma
                                  M0      = 4.5,
                                  b       = 1.0,
                                  Mmax    = NULL,
                                  T0,                             # tempo de origem da previsao (mesma escala do tempo do catalogo)
                                  horizonte,                      # duracao da janela de previsao [T0, T0+horizonte]
                                  mu_type = c("constant", "raster"),
                                  mu_const  = NULL,
                                  mu_raster = NULL,
                                  x_grid = NULL, y_grid = NULL,
                                  xlim = NULL, ylim = NULL,
                                  max_geracoes = 75,
                                  checar_criticidade = TRUE,
                                  semente = NULL
                                  ) {
    if (!is.null(semente)) set.seed(semente)
    mu_type <- match.arg(mu_type)
    if (checar_criticidade) verificar_criticidade(theta, b)
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
        warning("catalogo_historico contem eventos com time > T0; esses eventos foram descartados (nao sao 'historia conhecida' antes da previsao).")
        hist <- hist[hist$time <= T0, , drop = FALSE]
    }
    T_fim <- T0 + horizonte
    # ---- Amostradores auxiliares (autocontidos) ----
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
    # ---- 1. Réplicas futuras de eventos históricos ----
    # Para cada evento histórico i, o número de filhos que nascem em
    # (T0, T_fim) é Poisson com média kappa_i * [G(T_fim - t_i) - G(T0 - t_i)]
    # -- a diferença de CDFs garante que já não contamos os filhos que
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
    # ---- 2. Novos eventos de fundo na janela [T0, T_fim] ----
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
    # ---- 3. Junta as "sementes" (geração 0 da simulação futura) ----
    sementes <- rbind(bg_futuro, semente_hist)
    if (nrow(sementes) > 0) sementes$id <- seq_len(nrow(sementes))
    proximo_id <- nrow(sementes) + 1L
    # ---- 4. Propaga gerações a partir das sementes até T_fim ----
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
            warning("Numero maximo de geracoes atingido na simulacao de previsao -- verifique a razao de ramificacao antes de confiar no resultado.")
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
    # theta e' um "ensemble" (previsao deve amostrar theta a cada replica) se
    # for uma list() de vetores nomeados (ex.: thetas_ajustados_iran, vindo
    # do bootstrap). Um vetor nomeado unico como c(A=..., alpha=...) e'
    # atomico, nao list(), entao cai no ramo "theta unico" abaixo.
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
    # 'contagens' e' uma matriz n_bins x n_sim
    data.frame(
        bin_inicio = quebras[-length(quebras)],
        bin_fim    = quebras[-1],
        media      = rowMeans(contagens),
        mediana    = apply(contagens, 1, median),
        p025       = apply(contagens, 1, quantile, probs = 0.025),
        p975       = apply(contagens, 1, quantile, probs = 0.975)
    )
}

# ---- Agregação espacial: taxa esperada e P(>=1 evento) por célula do grid ----
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
        contagem_media    = contagem_media,   # numero esperado de eventos por celula
        prob_pelo_menos_1 = 1 - exp(-contagem_media)  # aproximacao de Poisson
    )
}

#====Predição Irã====

simulacoes_iran <- gerar_previsoes_etas(catalogo_historico = iran.cat,
                                        theta = thetas_iran[[length(thetas_iran)]],
                                        M0 = 4.5, b = 2, Mmax = 10,
                                        T0 = max(iran.cat$rtperiod), horizonte = 200,
                                        mu_type = "raster",
                                        mu_raster = rasters_iran[[length(rasters_iran)]],
                                        x_grid = x_grid_iran, y_grid = y_grid_iran,
                                        n_sim = 1000,
                                        semente = 123)

pred_temp_iran <- agregar_previsao_temporal(simulacoes_iran, 30)

pred_espac_iran <- agregar_previsao_espacial(simulacoes_iran, x_grid_iran, y_grid_iran)

pred_espac_iran$contagem_media %>% matrix(nrow = 128) %>% as.im() %>% plot()
