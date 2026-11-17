# =============================================================================
# Consistência assintótica do método (EM com suavização do background)
# Refeito a partir de estudo_consistencia.R com as correções acumuladas.
#
# Desenho: janela temporal crescente (n cresce com T), janela espacial fixa.
# Para cada cenário e cada nível de n, simula B catálogos INDEPENDENTES do mesmo
# theta verdadeiro, ajusta cada um a partir de um chute inicial fixo e resume
# viés, desvio-padrão, RMSE, cobertura de IC de Wald e recuperação do background.
#
# O que mudou em relação ao script antigo
#  1. Background verdadeiro paramétrico (exponenciais radiais + piso), em vez do
#     raster de outro ajuste (aj_peru$raster). Sem dependência de Inst_Peru.xlsx.
#  2. Funções atualizadas (patch_background.R): KDE por massa de célula com correção
#     de borda, mu nos eventos sem grade e com leave-one-out, raster inicial constante.
#  3. Catálogo simulado entra como data.frame (sem ETAS::catalog, POSIXct, nem
#     mag.threshold implícito); janela retangular aplicada explicitamente.
#  4. Seleção das réplicas por em_convergiu (EM atingiu o critério), não por
#     all(conv_codes == 0) && i < max_iter_em, que descartava réplicas só porque o
#     optim() falhou numa iteração intermediária (viés de seleção). A taxa de
#     convergência é reportada por nível.
#  5. Erros-padrão por dois métodos: Hessiana do passo M (subestima a incerteza)
#     e informação observada de Louis (informação completa menos a variância da
#     latente). Cobertura e razão SE/SD empírico para ambos.
#  6. Vazamento de borda medido em cada réplica (fração dos filhos de eventos
#     dentro da janela que caem fora). Com G = F = 1 na verossimilhança, esse
#     vazamento gera um viés assintótico em A que NÃO desaparece com T crescente
#     e janela fixa; o script o reporta e o compara ao viés observado.
#  7. Dois cenários: 'base' e 'borda_reduzida' (D menor => menos vazamento), para
#     separar viés de borda da consistência do restante do algoritmo.
#  8. Checagem de criticidade com o MESMO b da simulação.
#  9. Níveis definidos por n alvo (não por T), com paralelismo ajustado à memória
#     (cada processo guarda a lista de p_ij, n^2/2 valores) e retomada por checkpoint.
# 10. Métricas adicionais: erro L1 do background estimado e Brier de rho_j contra a
#     latente verdadeira.
#
# Requer functions.R e patch_background.R. Escrito sem execução em R no ambiente de
# redação: rode antes com NIVEIS_N pequeno, B_REPLICAS = 4 e NUM_CORES = 2.
# =============================================================================

suppressPackageStartupMessages({
    library(spatstat)   # nndist (d_j)
    library(parallel)
    library(dplyr)
})

#====Configurações====
ARQ_FUNCOES  <- "functions.R"
ARQ_PATCH    <- "patch_background.R"
ARQ_CHECK    <- "checkpoint_consistencia_novo.RData"
ARQ_SAIDA    <- "resultados_consistencia_novo.RData"
ARQ_PDF      <- "consistencia_assintotica.pdf"
RETOMAR      <- TRUE

M0           <- 4.5
B_GR         <- 2.0                       # valor-b da simulação (mesmo do script antigo)
MMAX         <- 10
XLIM         <- c(0, 8); YLIM <- c(0, 8)  # janela espacial fixa (graus)
N_GRID       <- 64                        # grade de análise (centros de célula)
N_GRID_SIM   <- 256                       # grade fina só para sortear o background
MU_TOT_DIA   <- 0.11                      # taxa total de background (eventos/dia na janela)
NIVEIS_N     <- c(500, 1000, 2000, 4000, 8000)   # n alvo (T é calculado a partir dele)
B_REPLICAS   <- c(100, 100, 60, 40, 20)          # réplicas por nível (reduza nos níveis grandes)
NUM_CORES    <- 15
RAM_GB       <- 64                        # memória disponível, para limitar o paralelismo
MULT_PARTIDA <- 1.4                       # chute inicial: theta0 = perturbação de 1.4 em todos os parâmetros
CALCULAR_SE  <- TRUE                      # Hessiana e Louis (exige guardar a lista de p_ij)
MAX_ITER_EM  <- 100

theta_base <- c("A" = 0.5, "alpha" = 1.5, "c" = 0.25,
                "p" = 1.5, "D" = 0.1, "q" = 2, "gamma" = 1)
theta_borda_reduzida <- replace(theta_base, "D", 0.01)
CENARIOS <- list(base = theta_base, borda_reduzida = theta_borda_reduzida)
nomes_theta <- names(theta_base)

lower_full <- c(A = 1e-6, alpha = 1e-6, c = 1e-6, p = 1.001, D = 1e-6, q = 1.001, gamma = 1e-6)
upper_full <- c(A = 15, alpha = 10, c = 15, p = 15, D = 15, q = 15, gamma = 10)

source(ARQ_FUNCOES)
source(ARQ_PATCH)
stopifnot(exists("simular_catalogo_etas"), exists("ajustar_etas_proprio"),
          exists("estimate_background_kernel_celula"), exists("log_lik_etas"),
          exists("gradiente_etas_completo"))

#====Background verdadeiro: exponenciais radiais + piso====
CENTROS    <- rbind(c(2.5, 5.5), c(5.5, 2.0))
AMPLITUDES <- c(1.0, 0.6)
ESCALAS    <- c(0.8, 1.2)
PISO_REL   <- 0.05
mu_forma <- function(x, y) {
    v <- rep(PISO_REL, length(x))
    for (k in seq_len(nrow(CENTROS))) {
        r <- sqrt((x - CENTROS[k, 1])^2 + (y - CENTROS[k, 2])^2)
        v <- v + AMPLITUDES[k] * exp(-r / ESCALAS[k])
    }
    v
}
centros_celula <- function(lim, N) lim[1] + (seq_len(N) - 0.5) * diff(lim) / N
raster_celula <- function(fun, lim_x, lim_y, Nx, Ny, sub = 5) {
    dx <- diff(lim_x) / Nx; dy <- diff(lim_y) / Ny
    off <- (seq_len(sub) - 0.5) / sub - 0.5
    xc <- centros_celula(lim_x, Nx); yc <- centros_celula(lim_y, Ny)
    Z <- matrix(0, Nx, Ny)
    for (a in off) for (b in off) Z <- Z + outer(xc + a * dx, yc + b * dy, fun)
    Z / sub^2
}
dxs <- diff(XLIM) / N_GRID_SIM; dys <- diff(YLIM) / N_GRID_SIM
I0  <- sum(raster_celula(mu_forma, XLIM, YLIM, N_GRID_SIM, N_GRID_SIM, sub = 2)) * dxs * dys
ESC <- MU_TOT_DIA / I0                              # integral de mu = MU_TOT_DIA eventos/dia
x_grid <- centros_celula(XLIM, N_GRID); y_grid <- centros_celula(YLIM, N_GRID)
raster_verd <- ESC * raster_celula(mu_forma, XLIM, YLIM, N_GRID, N_GRID)
x_sim <- centros_celula(XLIM, N_GRID_SIM); y_sim <- centros_celula(YLIM, N_GRID_SIM)
raster_sim <- ESC * raster_celula(mu_forma, XLIM, YLIM, N_GRID_SIM, N_GRID_SIM, sub = 2)

#====Funções auxiliares====
perturbar_theta <- function(theta, mult) {
    out <- theta
    for (nm in names(theta)) {
        out[nm] <- if (nm %in% c("p", "q")) 1 + (theta[nm] - 1) * mult else theta[nm] * mult
    }
    out
}

razao_ramificacao <- function(theta, b) {
    beta <- b * log(10)
    unname(theta["A"] * beta / (beta - theta["alpha"]))
}

# T que produz, em média, n_alvo eventos (regime estacionário: N ~ mu_tot * T / (1 - n))
T_para_n <- function(n_alvo, theta) n_alvo * (1 - razao_ramificacao(theta, B_GR)) / MU_TOT_DIA

# Nº de processos que cabem na memória (lista de p_ij: n^2/2 doubles, ~6 cópias em uso)
cores_para_n <- function(n) {
    max(1, min(NUM_CORES, floor(RAM_GB * 1e9 / (6 * 8 * n^2 / 2))))
}

# Vazamento de borda: fração dos filhos (já válidos no tempo) de eventos DENTRO da janela
# que caem FORA dela. Aproxima 1 - F_S médio ponderado pela produtividade.
na_janela <- function(x, y) x >= XLIM[1] & x <= XLIM[2] & y >= YLIM[1] & y <= YLIM[2]
vazamento_borda <- function(cat_full) {
    dentro <- na_janela(cat_full$x, cat_full$y)
    filhos <- which(cat_full$parent > 0)
    idp <- match(cat_full$parent[filhos], cat_full$id)
    pai_dentro <- dentro[idp]
    if (sum(pai_dentro) == 0) return(NA_real_)
    sum(pai_dentro & !dentro[filhos]) / sum(pai_dentro)
}

# Variância condicional do escore completo dada a latente (método de Louis).
# Condicional nos dados, u_j independentes com P(u_j = i) = p_ij e P(u_j = 0) = rho_j;
# o escore do termo de background em theta é zero.
var_escore_latente <- function(th, df, pij) {
    A <- unname(th["A"]); cc <- unname(th["c"]); pp <- unname(th["p"])
    D <- unname(th["D"]); qq <- unname(th["q"]); gam <- unname(th["gamma"])
    V <- matrix(0, 7, 7)
    for (j in 2:nrow(df)) {
        p <- pij[[j]]
        if (length(p) == 0) next
        i  <- 1:(j - 1)
        dt <- df$time[j] - df$time[i]
        d2 <- (df$x[j] - df$x[i])^2 + (df$y[j] - df$y[i])^2
        dm <- df$mag[i] - M0
        S_i <- D * exp(gam * dm); Vv <- d2 / S_i; U <- 1 + Vv; W <- qq * Vv / U - 1
        S <- cbind(rep(1 / A, length(i)), dm, (pp - 1) / cc - pp / (dt + cc),
                   1 / (pp - 1) + log(cc) - log(dt + cc), W / D, 1 / (qq - 1) - log(U), dm * W)
        m <- colSums(S * p)
        V <- V + crossprod(S * sqrt(p)) - tcrossprod(m)
    }
    dimnames(V) <- list(names(th), names(th))
    V
}

# Erros-padrão: (a) Hessiana do passo M; (b) informação observada de Louis
erros_padrao <- function(aj, df_ord) {
    th <- aj$theta
    Hc <- optimHess(th, log_lik_etas, gr = gradiente_etas_completo,
                    df_eventos = df_ord, matriz_p_ij = aj$probabilidades,
                    prob_total = aj$eventos$prob_disparo, raster_mu = aj$raster,
                    x_grid = x_grid, y_grid = y_grid, M0 = M0,
                    mu_eventos = aj$mu_eventos, integral_bg = sum(aj$rho_bg))
    dimnames(Hc) <- list(names(th), names(th))
    se_de <- function(I) {
        cv <- tryCatch(solve(I), error = function(e) NULL)
        if (is.null(cv) || any(!is.finite(diag(cv))) || any(diag(cv) <= 0)) return(rep(NA_real_, length(th)))
        setNames(sqrt(diag(cv)), names(th))
    }
    V <- var_escore_latente(th, df_ord, aj$probabilidades)
    list(se_hess = se_de(Hc), se_louis = se_de(Hc - V))
}

#====Uma réplica====
replica <- function(cen, idx_nivel, b) {
    theta_v <- CENARIOS[[cen]]
    T_niv   <- T_para_n(NIVEIS_N[idx_nivel], theta_v)
    semente <- 1e6 * match(cen, names(CENARIOS)) + 1e4 * idx_nivel + b
    t_ini <- Sys.time()
    tryCatch({
        cat_full <- simular_catalogo_etas(
            theta = theta_v, M0 = M0, b = B_GR, Mmax = MMAX, T_max = T_niv,
            mu_type = "raster", mu_raster = raster_sim, x_grid = x_sim, y_grid = y_sim,
            checar_criticidade = FALSE, semente = semente)
        leak <- vazamento_borda(cat_full)
        cat_w <- cat_full[na_janela(cat_full$x, cat_full$y), ]
        cat_w <- cat_w[order(cat_w$time), ]
        df <- data.frame(x = cat_w$x, y = cat_w$y, mag = cat_w$mag, time = cat_w$time)
        is_bg <- cat_w$parent == 0
        aj <- ajustar_etas_proprio(
            df, theta0 = perturbar_theta(theta_v, MULT_PARTIDA), M0 = M0,
            x_grid = x_grid, y_grid = y_grid, num_cores = 1, max_iter_em = MAX_ITER_EM,
            lower = lower_full, upper = upper_full, calcular_hessiana = FALSE,
            retornar_pij = CALCULAR_SE, verbose = FALSE)
        se <- if (CALCULAR_SE) erros_padrao(aj, aj$eventos[, c("x", "y", "mag", "time")]) else
            list(se_hess = rep(NA_real_, 7), se_louis = rep(NA_real_, 7))
        list(theta_hat = aj$theta, se_hess = se$se_hess, se_louis = se$se_louis,
             em_convergiu = aj$em_convergiu, todos_optim_ok = aj$todos_optim_ok,
             n_iter = aj$n_iter_em, n_eventos = nrow(df), T_max = T_niv, leak = leak,
             brier_rho = mean((aj$rho_bg - as.numeric(is_bg))^2),
             L1_bg = sum(abs(aj$raster - raster_verd)) / sum(raster_verd),
             tempo_min = as.numeric(difftime(Sys.time(), t_ini, units = "mins")), erro = NULL)
    }, error = function(e) list(erro = conditionMessage(e)))
}

#====Execução====
for (cen in names(CENARIOS)) verificar_criticidade(CENARIOS[[cen]], b = B_GR)
cat(sprintf("Razão de ramificação verdadeira: %.3f\n", razao_ramificacao(theta_base, B_GR)))

resultados <- setNames(lapply(names(CENARIOS), function(x) vector("list", length(NIVEIS_N))), names(CENARIOS))
if (RETOMAR && file.exists(ARQ_CHECK)) {
    load(ARQ_CHECK)
    cat("Retomando de", ARQ_CHECK, "\n")
}

for (cen in names(CENARIOS)) {
    for (l in seq_along(NIVEIS_N)) {
        if (!is.null(resultados[[cen]][[l]])) next
        nc <- cores_para_n(NIVEIS_N[l])
        cat(sprintf("\n==== Cenário %s | n alvo = %d (T = %.0f dias) | %d réplicas em %d processos ====\n",
                    cen, NIVEIS_N[l], T_para_n(NIVEIS_N[l], CENARIOS[[cen]]), B_REPLICAS[l], nc))
        ajustes <- mclapply(seq_len(B_REPLICAS[l]), function(b) replica(cen, l, b),
                            mc.cores = nc, mc.preschedule = FALSE)
        resultados[[cen]][[l]] <- list(n_alvo = NIVEIS_N[l], ajustes = ajustes)
        save(resultados, file = ARQ_CHECK, compress = TRUE)
    }
}

#====Análise====
falhou <- function(x) is.null(x) || inherits(x, "try-error") || !is.null(x$erro)
resumo <- do.call(rbind, lapply(names(CENARIOS), function(cen) {
    theta_v <- CENARIOS[[cen]]
    do.call(rbind, lapply(seq_along(NIVEIS_N), function(l) {
        aj_all <- resultados[[cen]][[l]]$ajustes
        sem_erro <- Filter(function(x) !falhou(x), aj_all)
        ok <- Filter(function(x) isTRUE(x$em_convergiu), sem_erro)
        if (length(ok) < 3) {
            warning(sprintf("Cenário %s, n = %d: menos de 3 réplicas válidas.", cen, NIVEIS_N[l]))
            return(NULL)
        }
        TH  <- do.call(rbind, lapply(ok, function(x) x$theta_hat)); colnames(TH) <- nomes_theta
        SEH <- do.call(rbind, lapply(ok, function(x) x$se_hess));   colnames(SEH) <- nomes_theta
        SEL <- do.call(rbind, lapply(ok, function(x) x$se_louis));  colnames(SEL) <- nomes_theta
        Bv  <- nrow(TH)
        vies <- colMeans(TH) - theta_v
        sdv  <- apply(TH, 2, sd)
        cob <- function(SE) sapply(nomes_theta, function(par) {
            ok_se <- is.finite(SE[, par])
            if (!any(ok_se)) return(NA_real_)
            mean(abs(TH[ok_se, par] - theta_v[par]) <= 1.96 * SE[ok_se, par])
        })
        razao <- function(SE) sapply(nomes_theta, function(par) mean(SE[, par], na.rm = TRUE) / sdv[par])
        data.frame(
            cenario = cen, nivel = l, n_alvo = NIVEIS_N[l],
            n_medio = mean(vapply(ok, function(x) x$n_eventos, numeric(1))),
            B_total = length(aj_all), B_valido = Bv,
            taxa_em_convergiu = length(ok) / length(sem_erro),
            frac_falhas = 1 - length(sem_erro) / length(aj_all),
            leak_medio = mean(vapply(ok, function(x) x$leak, numeric(1)), na.rm = TRUE),
            L1_bg_mediano = median(vapply(ok, function(x) x$L1_bg, numeric(1))),
            brier_mediano = median(vapply(ok, function(x) x$brier_rho, numeric(1))),
            parametro = nomes_theta, verdade = as.numeric(theta_v),
            vies = as.numeric(vies), vies_rel = as.numeric(vies / theta_v),
            sd = as.numeric(sdv), rmse = sqrt(as.numeric(vies)^2 + as.numeric(sdv)^2),
            z_vies = as.numeric(vies / (sdv / sqrt(Bv))),
            cob_hess = as.numeric(cob(SEH)), cob_louis = as.numeric(cob(SEL)),
            razao_se_hess = as.numeric(razao(SEH)), razao_se_louis = as.numeric(razao(SEL)),
            row.names = NULL)
    }))
}))

cat("\n===== Resumo por cenário, nível e parâmetro =====\n")
print(resumo, digits = 3, row.names = FALSE)

# Taxas de convergência: SD ~ n^(-1/2) (inclinação -0.5 em log-log) vale mesmo com viés persistente;
# MSE ~ 1/n só se o viés desaparecer. Reportar as duas separa variância de viés.
taxas <- do.call(rbind, lapply(split(resumo, list(resumo$cenario, resumo$parametro), drop = TRUE), function(d) {
    if (nrow(d) < 3) return(NULL)
    data.frame(cenario = d$cenario[1], parametro = d$parametro[1],
               inclinacao_sd = unname(coef(lm(log(sd) ~ log(n_medio), data = d))[2]),
               inclinacao_mse = unname(coef(lm(log(rmse^2) ~ log(n_medio), data = d))[2]),
               vies_rel_no_maior_n = d$vies_rel[which.max(d$n_medio)],
               z_vies_no_maior_n = d$z_vies[which.max(d$n_medio)])
}))
cat("\n===== Inclinações em log-log (esperado: SD ~ -0.5; MSE ~ -1 se viés -> 0) =====\n")
print(taxas, digits = 3, row.names = FALSE)

# Viés de A previsto pelo vazamento de borda (ordem de grandeza: -leak)
cat("\n===== Viés relativo de A observado vs. previsto pelo vazamento (-leak) =====\n")
print(subset(resumo, parametro == "A")[, c("cenario", "n_alvo", "vies_rel", "leak_medio", "z_vies")], digits = 3, row.names = FALSE)

#====Figuras====
pdf(ARQ_PDF, width = 11, height = 8)
cores_cob <- c(hess = "#B2182B", louis = "#2166AC")
for (cen in names(CENARIOS)) {
    dc <- subset(resumo, cenario == cen)
    if (nrow(dc) == 0) next
    # Viés relativo com IC de Monte Carlo
    par(mfrow = c(3, 3), mar = c(4, 4.2, 3, 1), oma = c(0, 0, 2, 0))
    for (par_nome in nomes_theta) {
        d <- subset(dc, parametro == par_nome)
        ep <- 1.96 * d$sd / sqrt(d$B_valido) / abs(d$verdade)
        plot(d$n_medio, d$vies_rel, type = "b", pch = 16, log = "x", xlab = "n médio de eventos",
             ylab = "Viés relativo", main = par_nome, ylim = range(c(d$vies_rel - ep, d$vies_rel + ep, 0)))
        arrows(d$n_medio, d$vies_rel - ep, d$n_medio, d$vies_rel + ep, angle = 90, code = 3, length = 0.03)
        abline(h = 0, lty = 2, col = "grey50")
        if (par_nome == "A") lines(d$n_medio, -d$leak_medio, col = 2, lty = 3, lwd = 2)
    }
    plot.new(); legend("center", c("Viés relativo (IC 95% de Monte Carlo)", "-vazamento de borda (previsão p/ A)"),
                       lty = c(1, 3), col = c(1, 2), bty = "n")
    mtext(sprintf("Cenário %s: viés relativo", cen), outer = TRUE, font = 2)

    # RMSE e SD em log-log
    par(mfrow = c(3, 3), mar = c(4, 4.2, 3, 1), oma = c(0, 0, 2, 0))
    for (par_nome in nomes_theta) {
        d <- subset(dc, parametro == par_nome)
        plot(d$n_medio, d$rmse, type = "b", pch = 16, log = "xy", xlab = "n médio de eventos",
             ylab = "RMSE (preto) e SD (azul)", main = par_nome, ylim = range(c(d$rmse, d$sd)))
        lines(d$n_medio, d$sd, type = "b", pch = 1, col = "#2166AC")
        lines(d$n_medio, d$sd[1] * (d$n_medio / d$n_medio[1])^(-0.5), lty = 2, col = "grey40")
    }
    plot.new(); legend("center", c("RMSE", "SD empírico", "referência n^(-1/2)"),
                       lty = c(1, 1, 2), pch = c(16, 1, NA), col = c(1, "#2166AC", "grey40"), bty = "n")
    mtext(sprintf("Cenário %s: RMSE e SD", cen), outer = TRUE, font = 2)

    # Cobertura e razão SE/SD
    par(mfrow = c(3, 3), mar = c(4, 4.2, 3, 1), oma = c(0, 0, 2, 0))
    for (par_nome in nomes_theta) {
        d <- subset(dc, parametro == par_nome)
        plot(d$n_medio, d$cob_hess, type = "b", pch = 16, col = cores_cob["hess"], log = "x", ylim = c(0, 1),
             xlab = "n médio de eventos", ylab = "Cobertura do IC 95%", main = par_nome)
        lines(d$n_medio, d$cob_louis, type = "b", pch = 17, col = cores_cob["louis"])
        abline(h = 0.95, lty = 2)
    }
    plot.new(); legend("center", c("Hessiana do passo M", "Louis", "nominal 0.95"),
                       lty = c(1, 1, 2), pch = c(16, 17, NA), col = c(cores_cob, "black"), bty = "n")
    mtext(sprintf("Cenário %s: cobertura de Wald", cen), outer = TRUE, font = 2)

    par(mfrow = c(3, 3), mar = c(4, 4.2, 3, 1), oma = c(0, 0, 2, 0))
    for (par_nome in nomes_theta) {
        d <- subset(dc, parametro == par_nome)
        plot(d$n_medio, d$razao_se_hess, type = "b", pch = 16, col = cores_cob["hess"], log = "x",
             ylim = range(c(d$razao_se_hess, d$razao_se_louis, 1), na.rm = TRUE),
             xlab = "n médio de eventos", ylab = "SE médio / SD empírico", main = par_nome)
        lines(d$n_medio, d$razao_se_louis, type = "b", pch = 17, col = cores_cob["louis"])
        abline(h = 1, lty = 2)
    }
    plot.new(); legend("center", c("Hessiana do passo M", "Louis", "calibrado = 1"),
                       lty = c(1, 1, 2), pch = c(16, 17, NA), col = c(cores_cob, "black"), bty = "n")
    mtext(sprintf("Cenário %s: calibração dos erros-padrão", cen), outer = TRUE, font = 2)
}

# Convergência, vazamento, background e latente por nível
par(mfrow = c(2, 2), mar = c(4, 4.2, 3, 1), oma = c(0, 0, 0, 0))
dn <- unique(resumo[, c("cenario", "n_medio", "taxa_em_convergiu", "leak_medio", "L1_bg_mediano", "brier_mediano")])
cc <- c(base = "black", borda_reduzida = "#D95F02")
for (campo in c("taxa_em_convergiu", "leak_medio", "L1_bg_mediano", "brier_mediano")) {
    plot(NA, xlim = range(dn$n_medio), ylim = range(dn[[campo]], na.rm = TRUE), log = "x",
         xlab = "n médio de eventos", ylab = campo, main = campo)
    for (cen in names(CENARIOS)) {
        d <- subset(dn, cenario == cen)
        if (nrow(d)) lines(d$n_medio, d[[campo]], type = "b", pch = 16, col = cc[cen])
    }
    legend("topright", names(cc), col = cc, pch = 16, bty = "n")
}
dev.off()

#====Salvar====
save(resultados, resumo, taxas, CENARIOS, NIVEIS_N, B_REPLICAS, MULT_PARTIDA, raster_verd,
     file = ARQ_SAIDA, compress = TRUE)
cat("\nResultados em", ARQ_SAIDA, "| figuras em", ARQ_PDF, "\n")
