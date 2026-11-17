# =============================================================================
# PATCH: background sem dependência da resolução da grade
#
# Problemas tratados
#  (1) d_min vs. espaçamento da grade: núcleos com h_j < Delta ficam mal
#      amostrados quando o KDE é avaliado apenas nos nós da grade.
#  (2) Hipótese  T * integral(mu_hat) = soma(rho_j): passa a valer EXATAMENTE
#      na janela (retângulo da grade), por construção.
#
# O que muda
#  - estimate_background_kernel_celula(): KDE por MASSA DE CÉLULA (diferenças de
#    pnorm nas bordas das células) + correção de borda c_j. Conserva a massa
#    para qualquer h_j, inclusive h_j << Delta. Vetorizado por produto de matrizes.
#  - mu_em_eventos(): mu_hat avaliado EXATAMENTE nos eventos (sem grade), com
#    opção leave-one-out (exclui a autocontribuição do evento j).
#  - calcular_probabilidades() e log_lik_etas(): aceitam o vetor mu_eventos
#    (e integral_bg), garantindo o MESMO mu nos passos E e M.
#  - ajustar_etas_proprio(): laço unificado usando as funções acima.
#
# Compatibilidade
#  - O raster devolvido continua na mesma grade e nas mesmas unidades
#    (eventos / unidade de tempo / unidade de área). Mapas, simulador e
#    resíduos continuam funcionando com ajuste$raster.
#  - Se mu_eventos = NULL, as funções usam func_nu() como antes.
#
# Convenção de grade: x_grid e y_grid são CENTROS de célula, com espaçamento
# constante; a janela é [x_1 - dx/2, x_N + dx/2] x [y_1 - dy/2, y_M + dy/2].
#
# Observação: este código foi escrito sem execução em R no ambiente de
# redação; a álgebra do estimador (conservação de massa, produto de matrizes)
# foi verificada numericamente em Python. Teste em um catálogo pequeno antes
# de rodar o catálogo completo.
# =============================================================================


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
                                 retornar_pij = FALSE,  # TRUE devolve a lista completa de p_ij (n^2/2 valores: pesada)
                                 calcular_hessiana = TRUE,  # FALSE acelera (útil em estudos de simulação/sensibilidade)
                                 verbose = TRUE
                                 ) {
    # `catalogo` pode ser um objeto de catálogo do ETAS ou um data.frame com colunas x, y, mag, time
    df_eventos <- if (is.data.frame(catalogo)) {
        catalogo[, c("x", "y", "mag", "time")]
    } else {
        extrair_eventos_catalog(catalogo, M0)
    }
    df_eventos <- df_eventos[order(df_eventos$time), ]
    rownames(df_eventos) <- NULL
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
                      control = list(maxit = 1000), hessian = calcular_hessiana)
        list(otim = otim, raster = bg$raster, mu_ev = mu_ev_novo)
    }

    thetas_hist  <- list(theta0)
    rasters_hist <- list(r0)
    theta_atual  <- theta0
    r_atual      <- r0
    mu_ev_atual  <- func_nu(r0, df_eventos, x_grid, y_grid)$lyr.1
    conv_codes   <- integer(0)
    hessiana_final <- NULL
    em_convergiu <- FALSE
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
        if (verbose) cat("====Iteração", i - 1, "finalizada====\n")
        if (convergiu_em(rasters_hist[[i]], rasters_hist[[i - 1]], tol_rel_raster, tol_abs_raster) &&
            convergiu_em(thetas_hist[[i]], thetas_hist[[i - 1]], tol_rel_theta, tol_abs_theta)) {
            em_convergiu <- TRUE
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
        em_convergiu       = em_convergiu,    # TRUE se os critérios de parada do EM foram atingidos (novo)
        codigo_optim_final = tail(conv_codes, 1),   # código do optim() no último passo M (0 = ok) (novo)
        todos_optim_ok     = all(conv_codes == 0),  # TRUE se TODOS os passos M terminaram com código 0 (novo)
        # convergiu: EM convergiu E o último passo M terminou bem. Falhas transitórias do optim em
        # iterações iniciais não reprovam o ajuste (ver todos_optim_ok).
        convergiu          = em_convergiu && tail(conv_codes, 1) == 0,
        thetas_hist        = if (guardar_historico) thetas_hist  else NULL,
        rasters_hist       = if (guardar_historico) rasters_hist else NULL
    )
}


# -----------------------------------------------------------------------------
# 6. Diagnóstico do background (para o texto e para a análise de sensibilidade)
# -----------------------------------------------------------------------------
diagnosticar_background <- function(ajuste, df_eventos, x_grid, y_grid, M0 = 4.5, num_cores = 1) {
    dx <- diff(x_grid)[1]; dy <- diff(y_grid)[1]
    h  <- ajuste$bw
    rho_bg <- ajuste$rho_bg          # já calculado no E-step final do ajuste
    mu_bilinear <- func_nu(ajuste$raster, df_eventos, x_grid, y_grid)$lyr.1
    list(
        espacamento           = c(dx = dx, dy = dy),
        frac_h_menor_meio_dx  = mean(h < max(dx, dy) / 2),   # núcleos mal resolvidos pela grade
        frac_h_no_piso        = mean(h <= min(h) + 1e-12),   # eventos no piso d_min
        razao_massa           = sum(ajuste$raster) * dx * dy / (sum(rho_bg) / ajuste$T_total),  # deve ser 1
        razao_mu_exato_bilin  = quantile(ajuste$mu_eventos / pmax(mu_bilinear, 1e-15),
                                         probs = c(0.05, 0.25, 0.5, 0.75, 0.95)),
        frac_background       = mean(rho_bg)
    )
}


# =============================================================================
# ALTERAÇÕES PONTUAIS nas funções de perfil (para usar o mesmo mu do ajuste)
# =============================================================================
# perfilar_rapido():
#   passo_e <- calcular_probabilidades(df_eventos, ajuste$theta, M0, ajuste$raster,
#                                      x_grid, y_grid, num_cores,
#                                      mu_eventos = ajuste$mu_eventos)
#   ... -log_lik_etas(theta_v, df_eventos, passo_e$probabilidades, passo_e$prob_total,
#                     ajuste$raster, x_grid, y_grid, M0 = M0,
#                     mu_eventos = ajuste$mu_eventos,
#                     integral_bg = sum(pmax(0, 1 - passo_e$prob_total)))
#
# perfilar_completo():
#   passo_e <- calcular_probabilidades(df_eventos, theta_ponto, M0, r_fixo,
#                                      x_grid, y_grid, num_cores_e,
#                                      mu_eventos = ajuste$mu_eventos)
#   obj <- function(theta_livre) { ...
#       log_lik_etas(theta_completo, df_eventos, passo_e$probabilidades, passo_e$prob_total,
#                    r_fixo, x_grid, y_grid, M0 = M0,
#                    mu_eventos = ajuste$mu_eventos,
#                    integral_bg = sum(pmax(0, 1 - passo_e$prob_total)))
#   }
#   (grad_obj não muda: gradiente_etas_completo não usa mu)
#
# Bootstrap paramétrico / resíduos: onde houver chamada a calcular_probabilidades()
# ou log_lik_etas() com ajuste$raster, passe também mu_eventos = ajuste$mu_eventos
# (ou, para catálogos simulados, o mu_eventos devolvido por ajustar_etas_proprio()).
# =============================================================================
