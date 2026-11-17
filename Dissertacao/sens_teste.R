# =============================================================================
# Estudo de sensibilidade ao chute inicial do EM (ETAS espaço-temporal)
# com background VERDADEIRO dado por uma função paramétrica (exponenciais radiais).
#
# Pergunta: partindo de chutes iniciais diferentes (theta0 e/ou background inicial),
# o EM converge ao MESMO ótimo? Se não, é ótimo local ou cume plano (identificabilidade)?
#
# Diferenças em relação a estudo_sensibilidade.R
#  - O background que gera os dados vem de uma função (soma de exponenciais radiais
#    mais um piso), não de um raster ajustado a outro catálogo. Isso permite também
#    medir o quanto o mu estimado recupera a verdade.
#  - Não roda o ajuste do catálogo do Irã nem usa ETAS::catalog(): o catálogo simulado
#    entra direto como data.frame (x, y, mag, time).
#  - Compara cada ajuste com o ajuste de REFERÊNCIA (partindo da verdade), e não só
#    com theta verdadeiro (que mistura erro amostral com sensibilidade numérica),
#    e usa a log-verossimilhança observada para separar ótimo local de cume plano.
#  - Inclui um terceiro fator: o chute inicial do background (constante em várias
#    escalas, e exponenciais radiais mal posicionadas).
#  - Paraleliza sobre os pontos de partida (processos externos), com o passo E
#    em 1 núcleo, o que escala melhor que paralelizar o passo E.
#  - Mantém ajustes cujo optim() falhou em iterações intermediárias, desde que o EM
#    tenha convergido (o critério antigo descartava esses casos: viés de seleção).
#
# Requer: functions.R (suas funções) e patch_background.R (versão atualizada).
# Escrito sem execução em R no ambiente de redação: rode primeiro com
# N_BG_ALVO pequeno, N_MULTIFATOR = 5 e MULT_OAT curto para conferir.
# =============================================================================

suppressPackageStartupMessages({
    library(spatstat)   # nndist (usado por d_j)
    library(parallel)
    library(dplyr)
})

#====Configurações====
ARQ_FUNCOES   <- "functions.R"
ARQ_PATCH     <- "patch_background.R"
ARQ_CHECK     <- "checkpoint_sensibilidade_param.RData"
ARQ_SAIDA     <- "resultados_sensibilidade_param.RData"
ARQ_PDF       <- "sensibilidade_param.pdf"

M0            <- 4.5
B_SIM         <- 1.5          # valor-b usado para simular magnitudes
MMAX_SIM      <- 10
T_SIM         <- 3650         # dias
XLIM          <- c(0, 8)      # janela retangular (graus)
YLIM          <- c(0, 8)
N_GRID        <- 64           # grade de análise (centros de célula)
N_GRID_SIM    <- 256          # grade fina usada só para sortear o background
N_BG_ALVO     <- 400          # nº esperado de eventos de background em T_SIM (controla o tamanho do catálogo)
NUM_CORES     <- 15           # processos externos (um ponto de partida por processo)
MAX_ITER_EM   <- 100
N_MULTIFATOR  <- 50
MULT_OAT      <- seq(0.25, 2, length.out = 10)
FAIXA_MULT    <- c(0.25, 4)
SEMENTE_SIM   <- 123
SEMENTE_LHS   <- 321
LOGLIK_TOL    <- 2           # diferença de log-verossimilhança considerada relevante
ERRO_TOL      <- 0.10        # distância relativa máxima ao ajuste de referência para "mesmo ótimo"

# Theta verdadeiro (mantido do script original)
theta_verdadeiro <- c("A" = 0.5, "alpha" = 1.5, "c" = 0.25,
                      "p" = 1.5, "D" = 0.1, "q" = 2, "gamma" = 1)
nomes_theta <- names(theta_verdadeiro)

source(ARQ_FUNCOES)
source(ARQ_PATCH)   # sobrescreve calcular_probabilidades, log_lik_etas, ajustar_etas_proprio etc.
stopifnot(exists("simular_catalogo_etas"), exists("ajustar_etas_proprio"),
          exists("estimate_background_kernel_celula"), exists("d_j"))

#====Background verdadeiro: soma de exponenciais radiais + piso====
CENTROS    <- rbind(c(2.5, 5.5), c(5.5, 2.0))   # centros (x, y)
AMPLITUDES <- c(1.0, 0.6)                       # amplitudes relativas
ESCALAS    <- c(0.8, 1.2)                       # escalas ell (graus): exp(-r / ell)
PISO_REL   <- 0.05                              # piso uniforme relativo ao pico

mu_forma <- function(x, y) {
    v <- rep(PISO_REL, length(x))
    for (k in seq_len(nrow(CENTROS))) {
        r <- sqrt((x - CENTROS[k, 1])^2 + (y - CENTROS[k, 2])^2)
        v <- v + AMPLITUDES[k] * exp(-r / ESCALAS[k])
    }
    v
}

centros_celula <- function(lim, N) lim[1] + (seq_len(N) - 0.5) * diff(lim) / N

# Média da função em cada célula (sub x sub pontos): evita sub-amostrar picos estreitos
raster_celula <- function(fun, lim_x, lim_y, Nx, Ny, sub = 5) {
    dx <- diff(lim_x) / Nx; dy <- diff(lim_y) / Ny
    off <- (seq_len(sub) - 0.5) / sub - 0.5
    xc <- centros_celula(lim_x, Nx); yc <- centros_celula(lim_y, Ny)
    Z <- matrix(0, Nx, Ny)
    for (a in off) for (b in off) Z <- Z + outer(xc + a * dx, yc + b * dy, fun)
    Z / sub^2
}

# Escala para que a integral de mu * T seja N_BG_ALVO
dxs <- diff(XLIM) / N_GRID_SIM; dys <- diff(YLIM) / N_GRID_SIM
I0  <- sum(raster_celula(mu_forma, XLIM, YLIM, N_GRID_SIM, N_GRID_SIM, sub = 2)) * dxs * dys
ESC <- (N_BG_ALVO / T_SIM) / I0
mu_verd <- function(x, y) ESC * mu_forma(x, y)    # eventos / dia / grau^2

x_grid <- centros_celula(XLIM, N_GRID); y_grid <- centros_celula(YLIM, N_GRID)
dx <- diff(x_grid)[1]; dy <- diff(y_grid)[1]
raster_verd <- ESC * raster_celula(mu_forma, XLIM, YLIM, N_GRID, N_GRID)      # verdade na grade de análise
x_sim <- centros_celula(XLIM, N_GRID_SIM); y_sim <- centros_celula(YLIM, N_GRID_SIM)
raster_sim <- ESC * raster_celula(mu_forma, XLIM, YLIM, N_GRID_SIM, N_GRID_SIM, sub = 2)
cat(sprintf("Background verdadeiro: integral*T = %.1f eventos (alvo %d)\n",
            sum(raster_verd) * dx * dy * T_SIM, N_BG_ALVO))

#====Simulação do catálogo (UM catálogo fixo)====
# Atenção: a checagem de criticidade usa o MESMO b da simulação (no script original
# usava b = 2.0 enquanto a simulação usava b = 1.5).
verificar_criticidade(theta_verdadeiro, b = B_SIM)

cat_sim <- simular_catalogo_etas(
    theta = theta_verdadeiro, M0 = M0, b = B_SIM, Mmax = MMAX_SIM, T_max = T_SIM,
    mu_type = "raster", mu_raster = raster_sim, x_grid = x_sim, y_grid = y_sim,
    checar_criticidade = FALSE, semente = SEMENTE_SIM
)
dentro  <- cat_sim$x >= XLIM[1] & cat_sim$x <= XLIM[2] & cat_sim$y >= YLIM[1] & cat_sim$y <= YLIM[2]
cat_sim <- cat_sim[dentro, ]
cat_sim <- cat_sim[order(cat_sim$time), ]
rownames(cat_sim) <- NULL
df_sim <- data.frame(x = cat_sim$x, y = cat_sim$y, mag = cat_sim$mag,
                     time = cat_sim$time - min(cat_sim$time))
is_bg_true <- cat_sim$parent == 0          # verdade da variável latente (disparados com pai fora da janela contam como disparados)
n_sim <- nrow(df_sim); T_obs <- max(df_sim$time)
cat(sprintf("Catálogo simulado: %d eventos | fração de background verdadeira = %.3f\n", n_sim, mean(is_bg_true)))
if (n_sim < 500 || n_sim > 8000)
    warning("Tamanho do catálogo fora de 500-8000: ajuste N_BG_ALVO (custo O(n^2) por iteração).")

d_j_sim <- d_j(df_sim[, c("x", "y")])

#====Log-verossimilhança observada (G = F = 1), para separar ótimo local de cume plano====
loglik_observada <- function(theta, df, mu_ev, raster, T_total) {
    lam <- vapply(seq_len(nrow(df)), function(j) {
        if (j == 1) return(mu_ev[j])
        i <- 1:(j - 1)
        mu_ev[j] + sum(kappa_func(df$mag[i], M0, theta["A"], theta["alpha"]) *
                       g_func(df$time[j] - df$time[i], theta["c"], theta["p"]) *
                       f_spatial((df$x[j] - df$x[i])^2 + (df$y[j] - df$y[i])^2, df$mag[i], M0,
                                 theta["D"], theta["q"], theta["gamma"]))
    }, numeric(1))
    sum(log(pmax(lam, 1e-300))) -
        (sum(raster) * dx * dy * T_total + sum(kappa_func(df$mag, M0, theta["A"], theta["alpha"])))
}
ll_oraculo <- loglik_observada(theta_verdadeiro, df_sim, mu_verd(df_sim$x, df_sim$y), raster_verd, T_obs)
cat(sprintf("Log-verossimilhança observada no oráculo (theta e mu verdadeiros): %.2f\n", ll_oraculo))

#====Design dos pontos de partida====
perturbar_theta <- function(theta, mult) {
    stopifnot(all(names(mult) == names(theta)))
    out <- theta
    for (nm in names(theta)) {
        out[nm] <- if (nm %in% c("p", "q")) 1 + (theta[nm] - 1) * mult[nm] else theta[nm] * mult[nm]
    }
    out
}
gerar_lhs <- function(n, k, semente = NULL) {
    if (!is.null(semente)) set.seed(semente)
    U <- matrix(runif(n * k), n, k)
    for (j in seq_len(k)) U[, j] <- (sample.int(n) - 1 + U[, j]) / n
    U
}

lista_theta <- list(); meta <- list()
adicionar <- function(theta0, tipo, par, mult, mu0) {
    lista_theta[[length(lista_theta) + 1]] <<- theta0
    meta[[length(meta) + 1]] <<- data.frame(tipo = tipo, parametro_variado = par, multiplicador = mult,
                                            mu0 = mu0, stringsAsFactors = FALSE)
}
# (a) referência: começa na verdade
adicionar(theta_verdadeiro, "REF", NA_character_, 1, "const_media")
# (b) OAT: um parâmetro por vez
for (par in nomes_theta) for (m in MULT_OAT) {
    mult <- setNames(rep(1, length(nomes_theta)), nomes_theta); mult[par] <- m
    adicionar(perturbar_theta(theta_verdadeiro, mult), "OAT", par, m, "const_media")
}
# (c) LHS multifator em escala log
U <- gerar_lhs(N_MULTIFATOR, length(nomes_theta), SEMENTE_LHS); colnames(U) <- nomes_theta
M_lhs <- exp(log(FAIXA_MULT[1]) + U * (log(FAIXA_MULT[2]) - log(FAIXA_MULT[1])))
for (r in seq_len(N_MULTIFATOR)) {
    adicionar(perturbar_theta(theta_verdadeiro, setNames(M_lhs[r, ], nomes_theta)),
              "LHS", NA_character_, NA_real_, "const_media")
}
# (d) chute inicial do background (theta0 = verdade)
for (tp in c("const_x10", "const_x0.1", "uniforme_1", "exp_radial_centro", "exp_radial_canto")) {
    adicionar(theta_verdadeiro, "MU0", NA_character_, NA_real_, tp)
}
theta0_mat <- do.call(rbind, lista_theta); colnames(theta0_mat) <- nomes_theta
meta <- do.call(rbind, meta)
n_pontos <- nrow(meta)
cat("Pontos de partida:", n_pontos, "\n"); print(table(meta$tipo))

# Chutes iniciais do background. Todos normalizados para integral = n/T (todos os
# eventos como background), exceto as variantes de escala e o 'uniforme_1' (1 evento/dia/grau^2)
raster_inicial <- function(tipo) {
    Nx <- length(x_grid); Ny <- length(y_grid)
    media <- n_sim / (T_obs * Nx * dx * Ny * dy)
    forma <- function(fun) { Z <- outer(x_grid, y_grid, fun); Z * (n_sim / T_obs) / (sum(Z) * dx * dy) }
    switch(tipo,
           const_media       = matrix(media, Nx, Ny),
           const_x10         = matrix(10 * media, Nx, Ny),
           "const_x0.1"      = matrix(0.1 * media, Nx, Ny),
           uniforme_1        = matrix(1, Nx, Ny),
           exp_radial_centro = forma(function(x, y) exp(-sqrt((x - mean(XLIM))^2 + (y - mean(YLIM))^2) / 1.5)),
           exp_radial_canto  = forma(function(x, y) exp(-sqrt((x - XLIM[1])^2 + (y - YLIM[1])^2) / 1.5)),
           stop("tipo de mu0 desconhecido: ", tipo))
}

#====Ajuste de um ponto de partida====
ajustar_ponto <- function(k) {
    t0 <- theta0_mat[k, ]
    t_ini <- Sys.time()
    fit <- tryCatch(
        ajustar_etas_proprio(df_sim, t0, M0, x_grid, y_grid, d_j_pre = d_j_sim,
                             r0 = raster_inicial(meta$mu0[k]), num_cores = 1,
                             max_iter_em = MAX_ITER_EM, calcular_hessiana = FALSE, verbose = FALSE),
        error = function(e) list(erro = conditionMessage(e)))
    if (!is.null(fit$erro)) return(list(erro = fit$erro))
    ll <- loglik_observada(fit$theta, df_sim, fit$mu_eventos, fit$raster, fit$T_total)
    list(theta = fit$theta, raster = fit$raster, loglik = ll,
         em_convergiu = fit$em_convergiu, codigo_optim_final = fit$codigo_optim_final,
         todos_optim_ok = fit$todos_optim_ok, n_iter = fit$n_iter_em,
         brier_rho = mean((fit$rho_bg - as.numeric(is_bg_true))^2),
         L1_bg = sum(abs(fit$raster - raster_verd)) / sum(raster_verd),
         corr_bg = cor(as.vector(fit$raster), as.vector(raster_verd)),
         tempo_min = as.numeric(difftime(Sys.time(), t_ini, units = "mins")))
}

#====Execução (lotes, com checkpoint)====
resultados <- vector("list", n_pontos)
lotes <- split(seq_len(n_pontos), ceiling(seq_len(n_pontos) / (2 * NUM_CORES)))
for (b in seq_along(lotes)) {
    idx <- lotes[[b]]
    res <- mclapply(idx, ajustar_ponto, mc.cores = NUM_CORES, mc.preschedule = FALSE)
    resultados[idx] <- res
    save(resultados, meta, theta0_mat, theta_verdadeiro, file = ARQ_CHECK, compress = TRUE)
    cat(sprintf("Lote %d/%d concluído (%d/%d pontos)\n", b, length(lotes), max(idx), n_pontos))
}

#====Análise====
falhou <- function(r) is.null(r) || inherits(r, "try-error") || !is.null(r$erro)
ref <- resultados[[1]]
if (falhou(ref)) stop("O ajuste de referência (partindo da verdade) falhou; sem ele não há comparação.")

linhas <- lapply(seq_len(n_pontos), function(k) {
    r <- resultados[[k]]
    if (falhou(r)) return(data.frame(k = k, falha = TRUE))
    rel <- (r$theta - ref$theta) / abs(ref$theta)
    data.frame(k = k, falha = FALSE, em_convergiu = r$em_convergiu, optim_final_ok = r$codigo_optim_final == 0,
               todos_optim_ok = r$todos_optim_ok, n_iter = r$n_iter, loglik = r$loglik,
               dloglik = r$loglik - ref$loglik,
               dist_ref_max = max(abs(rel)), dist_ref_media = mean(abs(rel)),
               brier_rho = r$brier_rho, L1_bg = r$L1_bg, corr_bg = r$corr_bg, tempo_min = r$tempo_min)
})
tab <- do.call(rbind, lapply(linhas, function(d) {
    cols <- c("k", "falha", "em_convergiu", "optim_final_ok", "todos_optim_ok", "n_iter", "loglik", "dloglik",
              "dist_ref_max", "dist_ref_media", "brier_rho", "L1_bg", "corr_bg", "tempo_min")
    for (cc in setdiff(cols, names(d))) d[[cc]] <- NA
    d[, cols]
}))
tab <- cbind(meta, tab[, -1])
tab$status <- with(tab, ifelse(falha, "falha",
                        ifelse(!em_convergiu, "EM_nao_convergiu",
                        ifelse(dloglik < -LOGLIK_TOL, "otimo_local",
                        ifelse(dloglik > LOGLIK_TOL, "melhor_que_ref",
                        ifelse(dist_ref_max > ERRO_TOL, "cume_plano", "mesmo_otimo"))))))

cat("\n===== Referência (começa na verdade) vs. verdade =====\n")
print(round(rbind(verdade = theta_verdadeiro, ref = ref$theta, erro_rel = (ref$theta - theta_verdadeiro) / theta_verdadeiro), 4))
cat(sprintf("loglik: referência = %.2f | oráculo = %.2f\n", ref$loglik, ll_oraculo))

cat("\n===== Status dos ajustes por tipo de ponto de partida =====\n")
print(table(tab$tipo, tab$status))
cat("\n===== OAT: por parâmetro perturbado =====\n")
oat <- subset(tab, tipo == "OAT")
print(do.call(rbind, lapply(split(oat, oat$parametro_variado), function(d)
    data.frame(parametro = d$parametro_variado[1], frac_mesmo_otimo = mean(d$status == "mesmo_otimo"),
               dist_ref_mediana = median(d$dist_ref_max, na.rm = TRUE), dist_ref_max = max(d$dist_ref_max, na.rm = TRUE),
               n_iter_mediana = median(d$n_iter, na.rm = TRUE)))), digits = 3)
cat("\n===== LHS multifator =====\n")
lhs <- subset(tab, tipo == "LHS")
cat(sprintf("mesmo ótimo: %.0f%% | cume plano: %.0f%% | ótimo local: %.0f%% | EM não convergiu: %.0f%% | falha: %.0f%%\n",
            100 * mean(lhs$status == "mesmo_otimo"), 100 * mean(lhs$status == "cume_plano"),
            100 * mean(lhs$status == "otimo_local"), 100 * mean(lhs$status == "EM_nao_convergiu"),
            100 * mean(lhs$status == "falha")))
cat("\n===== Chute inicial do background (theta0 = verdade) =====\n")
print(subset(tab, tipo == "MU0")[, c("mu0", "status", "dist_ref_max", "dloglik", "n_iter", "L1_bg", "corr_bg")], digits = 3)
cat("\n===== Recuperação do background e da latente (todos os ajustes que convergiram) =====\n")
ok <- subset(tab, !falha & em_convergiu)
cat(sprintf("erro L1 relativo do mu estimado: mediana %.3f (min %.3f, max %.3f) | Brier de rho_j: mediana %.4f\n",
            median(ok$L1_bg), min(ok$L1_bg), max(ok$L1_bg), median(ok$brier_rho)))
cat(sprintf("Ajustes com optim() falhando em alguma iteração, mas EM convergido (mantidos): %d\n",
            sum(!tab$falha & tab$em_convergiu & !tab$todos_optim_ok, na.rm = TRUE)))

#====Figuras====
pdf(ARQ_PDF, width = 10, height = 7)
cores <- c(REF = "black", OAT = "#2166AC", LHS = "#1B9E77", MU0 = "#D95F02")
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
image(x_grid, y_grid, raster_verd, col = hcl.colors(100, "YlOrRd", rev = TRUE), asp = 1,
      xlab = "x", ylab = "y")
zr <- range(c(ref$raster, raster_verd))
image(x_grid, y_grid, ref$raster, col = hcl.colors(100, "YlOrRd", rev = TRUE), zlim = zr, asp = 1,
      xlab = "x", ylab = "y", main = "Background estimado (referência)")

par(mfrow = c(1, 2), mar = c(5, 4.5, 3, 1))
boxplot(dist_ref_max ~ tipo, data = subset(tab, !falha & tipo != "REF"), log = "y", col = "grey85", xlab = "",
        ylab = "Distância relativa máxima ao ajuste de referência", main = "Distância ao ótimo de referência")
abline(h = ERRO_TOL, lty = 2, col = 2)
plot(tab$dist_ref_max, tab$dloglik, col = cores[tab$tipo], pch = 16, log = "x",
     xlab = "Maximum relative distance to reference", ylab = "loglik - loglik(ref)")
abline(h = c(-LOGLIK_TOL, LOGLIK_TOL), lty = 2, col = 2); abline(v = ERRO_TOL, lty = 2, col = 2)
legend("topleft", names(cores), col = cores, pch = 16, bty = "n")

par(mfrow = c(3, 3), mar = c(4, 4.2, 2.5, 1))
for (par_nome in nomes_theta) {
    d <- subset(tab, tipo == "OAT" & parametro_variado == par_nome & !falha)
    plot(d$multiplicador, d$dist_ref_max, type = "b", pch = 16, log = "y", xlab = "Multiplicador do chute",
         ylab = "Dist. ao ajuste de referência", main = par_nome)
    abline(h = ERRO_TOL, lty = 2, col = 2)
}
dev.off()

#====Salvar====
save(theta_verdadeiro, theta0_mat, meta, resultados, tab, ll_oraculo, df_sim, is_bg_true,
     raster_verd, x_grid, y_grid, file = ARQ_SAIDA, compress = TRUE)
cat("\nResultados em", ARQ_SAIDA, "| figuras em", ARQ_PDF, "\n")
