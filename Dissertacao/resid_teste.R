# =============================================================================
# Bondade de ajuste e análise de resíduos do ETAS espaço-temporal (Peru)
# a partir de ajuste_peru.RData (objeto aj_peru).
#
# Entrada : ajuste_peru.RData  (aj_peru, t1_peru, t2_peru)
#           Inst_Peru.xlsx     (opcional, só para reconstruir a grade; ver abaixo)
# Saída   : residuos_peru.pdf  (figuras)  e  residuos_peru.RData  (res_peru)
#
# O que é feito
#  1. Compensador do processo marginal temporal e teorema da mudança aleatória
#     de tempo: tau_j, U_j, teste KS, Ljung-Box, gráficos.
#  2. Resíduos de Pearson por janela temporal.
#  3. Resíduos espaciais por pixel: contagens observadas vs. esperadas, com o
#     termo de disparo calculado por quadratura estratificada do kernel de
#     Pareto (trata o vazamento pela borda da janela), resíduo suavizado e
#     perfis marginais.
#  4. Diagnóstico de efeito de borda: rho_j vs. distância à borda.
#
# Pontos importantes
#  - O .RData NÃO guarda x_grid, y_grid nem as funções. A grade é reconstruída
#    com o mesmo código do ajuste (precisa de Inst_Peru.xlsx e do pacote ETAS);
#    sem o xlsx, usa-se a caixa dos eventos como aproximação (com aviso).
#  - O compensador usa a integral EXATA da intensidade na janela (fração do
#    kernel que cai dentro, F_S, por quadratura), e não a aproximação F_S = 1
#    usada na verossimilhança do ajuste. É o que o teorema da mudança de tempo
#    exige. A versão com F_S = 1 é calculada para comparação.
#  - Os p-valores do KS/Ljung-Box são aproximados: theta e mu foram estimados
#    dos mesmos dados (plug-in). Para p-valores formais, calibre por bootstrap
#    paramétrico.
#  - Escrito sem execução em R no ambiente de redação. Teste primeiro com
#    NUM_CORES = 1 e K_RAD/L_ANG menores para conferir.
# =============================================================================

suppressPackageStartupMessages({
    library(parallel)
    library(Matrix)      # sparseMatrix (acumulação rápida em pixels)
})

#====Configurações====
ARQ_RDATA    <- "~/ajuste_peru.RData"
ARQ_XLSX     <- "Inst_Peru.xlsx"     # NULL ou inexistente => grade aproximada
ARQ_PDF      <- "residuos_peru.pdf"
ARQ_SAIDA    <- "residuos_peru.RData"
M0           <- 4.5
NUM_CORES    <- max(1, parallel::detectCores() - 2)
N_BINS_TEMPO <- 60       # janelas temporais para o resíduo de Pearson
AGREGAR      <- 4        # fator de agregação dos pixels (64 -> 16)
E_MIN        <- 5        # esperado mínimo por pixel agregado
SIGMA_SUAV   <- 2        # largura (em células) do suavizador do resíduo
K_RAD        <- 16       # estratos radiais da quadratura do kernel
L_ANG        <- 24       # estratos angulares da quadratura do kernel
BLOCO        <- 400      # eventos por bloco na quadratura (controla memória)
COMPARAR_SEM_FS <- TRUE  # calcula também a versão com F_S = 1
SEMENTE      <- 2026

RNGkind("L'Ecuyer-CMRG")
set.seed(SEMENTE)

#====Funções do modelo====
kappa_func <- function(m, m0, A, alpha) A * exp(alpha * (m - m0))

# CDF do kernel temporal de Omori-Utsu (integral de g de 0 a s); zero para s <= 0
G_omori <- function(s, c, p) {
    out <- numeric(length(s))
    pos <- s > 0
    out[pos] <- 1 - (c / (s[pos] + c))^(p - 1)
    out
}

#====Carregar o ajuste====
load(ARQ_RDATA)
stopifnot(exists("aj_peru"))
campos <- c("theta", "raster", "eventos", "rho_bg", "T_total")
faltam <- setdiff(campos, names(aj_peru))
if (length(faltam) > 0) {
    stop("aj_peru não contém: ", paste(faltam, collapse = ", "),
         ". Reajuste com a versão de ajustar_etas_proprio que devolve rho_bg e eventos.")
}

th      <- aj_peru$theta
A_par   <- unname(th["A"]);  alpha_par <- unname(th["alpha"])
c_par   <- unname(th["c"]);  p_par     <- unname(th["p"])
D_par   <- unname(th["D"]);  q_par     <- unname(th["q"]);  gamma_par <- unname(th["gamma"])
raster  <- aj_peru$raster
ev      <- aj_peru$eventos
T_total <- aj_peru$T_total
n       <- nrow(ev)
Nx <- nrow(raster); Ny <- ncol(raster)

stopifnot(!is.unsorted(ev$time))   # eventos em ordem cronológica
cat(sprintf("Eventos: %d | T = %.1f dias (%.1f anos) | grade %d x %d\n",
            n, T_total, T_total / 365.25, Nx, Ny))
cat("Parâmetros estimados:\n"); print(round(th, 5))

#====Grade e janela====
reconstruir_grade <- function() {
    usa_xlsx <- !is.null(ARQ_XLSX) && file.exists(ARQ_XLSX) &&
        all(vapply(c("readxl", "dplyr", "ETAS"), requireNamespace, logical(1), quietly = TRUE))
    if (usa_xlsx) {
        suppressPackageStartupMessages({library(dplyr); library(readxl); library(ETAS)})
        peru <- read_xlsx(ARQ_XLSX) %>%
            select(date = `fecha UTC`, time = `hora UTC`,
                   lat = `latitud (º)`, long = `longitud (º)`,
                   prof = `profundidad (km)`, mag = `magnitud (M)`) %>%
            mutate(date = as.Date(date), lat = as.numeric(lat), long = as.numeric(long),
                   prof = as.numeric(prof), mag = as.numeric(mag)) %>%
            filter(mag >= M0)
        cat_peru <- catalog(peru)
        rl <- range(cat_peru$region.poly$long); rt <- range(cat_peru$region.poly$lat)
        n_cat <- sum(cat_peru$revents[, "flag"] %in% c(0, 1))
        if (n_cat != n) warning(sprintf("Catálogo reconstruído tem %d eventos; o ajuste usou %d.", n_cat, n))
        cat("Grade reconstruída a partir de", ARQ_XLSX, "\n")
    } else {
        rl <- range(ev$x); rt <- range(ev$y)
        warning("Inst_Peru.xlsx indisponível: grade aproximada pela caixa dos eventos. ",
                "Confirme que coincide com seq(range(region.poly), length.out = n) do ajuste.")
    }
    list(x = seq(rl[1], rl[2], length.out = Nx), y = seq(rt[1], rt[2], length.out = Ny))
}
gr <- reconstruir_grade()
x_grid <- gr$x; y_grid <- gr$y
dx <- diff(x_grid)[1]; dy <- diff(y_grid)[1]
# Convenção do ajuste: x_grid/y_grid são CENTROS de célula
bx <- c(x_grid - dx / 2, x_grid[Nx] + dx / 2)
by <- c(y_grid - dy / 2, y_grid[Ny] + dy / 2)
janela <- list(x0 = bx[1], x1 = bx[Nx + 1], y0 = by[1], y1 = by[Ny + 1])
stopifnot(all(ev$x >= janela$x0 & ev$x <= janela$x1 & ev$y >= janela$y0 & ev$y <= janela$y1))

# Massa do background na janela (identidade T * integral(mu) = soma(rho_j) no ajuste)
mu_tot <- sum(raster) * dx * dy
cat(sprintf("T * integral(mu) = %.2f | soma(rho_bg) = %.2f (devem coincidir)\n",
            mu_tot * T_total, sum(ev$rho_bg)))

#====Quadratura do kernel de Pareto: F_S e esperado espacial do disparo====
# Para cada evento i, o deslocamento dos filhos tem raio r com
#   r^2 = S_i [ (1-u)^(-1/(q-1)) - 1 ],  u ~ U(0,1),  e ângulo uniforme.
# Amostragem estratificada com jitter (K_RAD x L_ANG estratos, cada um com massa
# 1/(K*L)) dá a fração do kernel que cai em cada pixel e dentro da janela (F_S).
S_i <- D_par * exp(gamma_par * (ev$mag - M0))
w_i <- kappa_func(ev$mag, M0, A_par, alpha_par) * G_omori(T_total - ev$time, c_par, p_par)

quadratura_disparo <- function() {
    kk <- rep(0:(K_RAD - 1), each = L_ANG)
    ll <- rep(0:(L_ANG - 1), times = K_RAD)
    nq <- K_RAD * L_ANG
    blocos <- split(seq_len(n), ceiling(seq_len(n) / BLOCO))
    res <- mclapply(blocos, function(idx) {
        m  <- length(idx)
        U  <- (matrix(kk, m, nq, byrow = TRUE) + matrix(runif(m * nq), m, nq)) / K_RAD
        Th <- 2 * pi * (matrix(ll, m, nq, byrow = TRUE) + matrix(runif(m * nq), m, nq)) / L_ANG
        R  <- sqrt(S_i[idx] * ((1 - U)^(-1 / (q_par - 1)) - 1))
        X  <- ev$x[idx] + R * cos(Th)
        Y  <- ev$y[idx] + R * sin(Th)
        ix <- findInterval(X, bx)
        iy <- findInterval(Y, by)
        dentro <- matrix(ix >= 1 & ix <= Nx & iy >= 1 & iy <= Ny, m, nq)
        W  <- matrix(w_i[idx] / nq, m, nq)
        Mx <- as.matrix(sparseMatrix(i = ix[dentro], j = iy[dentro], x = W[dentro],
                                     dims = c(Nx, Ny)))
        list(M = Mx, FS = rowMeans(dentro))
    }, mc.cores = NUM_CORES, mc.set.seed = TRUE)
    list(M  = Reduce(`+`, lapply(res, `[[`, "M")),
         FS = unlist(lapply(res, `[[`, "FS"), use.names = FALSE))
}
cat("Quadratura do kernel espacial... ")
quad <- quadratura_disparo()
FS <- quad$FS
cat("ok\n")
cat(sprintf("F_S (fração do kernel dentro da janela): média %.3f | mínimo %.3f | eventos com F_S < 0.9: %.1f%%\n",
            mean(FS), min(FS), 100 * mean(FS < 0.9)))

#====1. Resíduos temporais (mudança aleatória de tempo)====
# Lambda(t) = mu_tot * t + sum_{t_i < t} kappa_i * F_S(i) * G(t - t_i)
Lambda_hat <- function(tt, kapF) {
    unlist(mclapply(tt, function(t) {
        k <- findInterval(t, ev$time, left.open = TRUE)   # nº de eventos com time < t
        if (k == 0) return(mu_tot * t)
        mu_tot * t + sum(kapF[1:k] * G_omori(t - ev$time[1:k], c_par, p_par))
    }, mc.cores = NUM_CORES), use.names = FALSE)
}
kap_i  <- kappa_func(ev$mag, M0, A_par, alpha_par)

diagnostico_tempo <- function(kapF, rotulo) {
    cat("Compensador nos tempos dos eventos (", rotulo, ")... ", sep = ""); flush.console()
    tau  <- Lambda_hat(ev$time, kapF)
    dtau <- diff(c(0, tau))
    if (any(dtau <= 0)) warning("dtau <= 0 em ", sum(dtau <= 0), " eventos (empates de tempo?).")
    U    <- 1 - exp(-pmax(dtau, 0))
    ks   <- suppressWarnings(ks.test(U, "punif"))
    lb   <- Box.test(U, lag = 20, type = "Ljung-Box")
    cat("ok\n")
    list(tau = tau, dtau = dtau, U = U, ks = ks, lb = lb)
}
tempo_fs <- diagnostico_tempo(kap_i * FS, "com F_S")
if (COMPARAR_SEM_FS) tempo_sem <- diagnostico_tempo(kap_i, "F_S = 1")

# Contagem total esperada vs observada (Lambda(T) deve ser próximo de n)
Lambda_T <- Lambda_hat(T_total, kap_i * FS)
z_total  <- (n - Lambda_T) / sqrt(Lambda_T)

# Resíduo de Pearson por janela temporal
bordas_t <- seq(0, T_total, length.out = N_BINS_TEMPO + 1)
Lam_bord <- Lambda_hat(bordas_t, kap_i * FS)
E_t <- diff(Lam_bord)
O_t <- hist(ev$time, breaks = bordas_t, plot = FALSE)$counts
R_t <- (O_t - E_t) / sqrt(E_t)
anos_centro <- (bordas_t[-1] + bordas_t[-length(bordas_t)]) / 2 / 365.25

#====2. Resíduos espaciais (pixels)====
ix_ev <- findInterval(ev$x, bx, all.inside = TRUE)
iy_ev <- findInterval(ev$y, by, all.inside = TRUE)
O_pix <- matrix(tabulate(ix_ev + (iy_ev - 1) * Nx, nbins = Nx * Ny), Nx, Ny)
E_bg  <- T_total * raster * dx * dy            # background: massa por célula
E_pix <- E_bg + quad$M                          # + disparo (quadratura)
cat(sprintf("Soma observada = %d | soma esperada nos pixels = %.1f (Lambda(T) = %.1f)\n",
            sum(O_pix), sum(E_pix), Lambda_T))

agregar <- function(M, f) {
    nx <- floor(nrow(M) / f) * f; ny <- floor(ncol(M) / f) * f
    M <- M[1:nx, 1:ny, drop = FALSE]
    out <- matrix(0, nx / f, ny / f)
    for (a in 1:f) for (b in 1:f) out <- out + M[seq(a, nx, by = f), seq(b, ny, by = f)]
    out
}
O_ag <- agregar(O_pix, AGREGAR); E_ag <- agregar(E_pix, AGREGAR)
x_ag <- colMeans(matrix(x_grid[1:(nrow(O_ag) * AGREGAR)], nrow = AGREGAR))
y_ag <- colMeans(matrix(y_grid[1:(ncol(O_ag) * AGREGAR)], nrow = AGREGAR))
R_pix <- ifelse(E_ag >= E_MIN, (O_ag - E_ag) / sqrt(E_ag), NA)
chi2_esp <- sum(R_pix^2, na.rm = TRUE); gl_esp <- sum(!is.na(R_pix))

# Suavização gaussiana separável (borda replicada)
blur2d <- function(M, sigma) {
    r <- ceiling(4 * sigma); k <- dnorm(-r:r, sd = sigma); k <- k / sum(k)
    conv1 <- function(v) {
        vv <- c(rep(v[1], r), v, rep(v[length(v)], r))
        as.numeric(stats::filter(vv, k, sides = 2))[(r + 1):(r + length(v))]
    }
    M1 <- t(apply(M, 1, conv1))   # ao longo de y
    apply(M1, 2, conv1)           # ao longo de x
}
res_suav <- blur2d(O_pix - E_pix, SIGMA_SUAV)
O_s <- blur2d(O_pix, SIGMA_SUAV); E_s <- blur2d(E_pix, SIGMA_SUAV)
razao_log <- ifelse(E_s >= 0.5, log((O_s + 0.5) / (E_s + 0.5)), NA)

#====3. Efeito de borda: rho_j vs distância à borda====
d_borda <- pmin(ev$x - janela$x0, janela$x1 - ev$x, ev$y - janela$y0, janela$y1 - ev$y)
cortes  <- c(0, 0.25, 0.5, 1, 2, 4, Inf)
faixa   <- cut(d_borda, cortes, include.lowest = TRUE)
rho_faixa <- tapply(ev$rho_bg, faixa, mean)
n_faixa   <- as.integer(table(faixa))
soma_rho_pix <- as.matrix(sparseMatrix(i = ix_ev, j = iy_ev, x = ev$rho_bg, dims = c(Nx, Ny)))
rho_ag <- ifelse(O_ag > 0, agregar(soma_rho_pix, AGREGAR) / O_ag, NA)

#====Resumo====
resumo <- list(
    n_eventos            = n,
    Lambda_T_vs_n        = c(Lambda_T = Lambda_T, n = n, z = z_total),
    KS_com_FS            = c(D = unname(tempo_fs$ks$statistic), p = tempo_fs$ks$p.value),
    LjungBox_com_FS      = c(X2 = unname(tempo_fs$lb$statistic), p = tempo_fs$lb$p.value),
    KS_F_S_igual_1       = if (COMPARAR_SEM_FS) c(D = unname(tempo_sem$ks$statistic), p = tempo_sem$ks$p.value) else NULL,
    Pearson_tempo        = c(media_R2 = mean(R_t^2), frac_abs_gt2 = mean(abs(R_t) > 2)),
    Pearson_espacial     = c(chi2 = chi2_esp, pixels = gl_esp, media_R2 = chi2_esp / gl_esp,
                             frac_abs_gt2 = mean(abs(R_pix) > 2, na.rm = TRUE)),
    FS                   = c(media = mean(FS), frac_menor_0.9 = mean(FS < 0.9)),
    rho_por_dist_borda   = data.frame(faixa = levels(faixa), n = n_faixa, rho_medio = as.numeric(rho_faixa)),
    frac_background      = mean(ev$rho_bg)
)
cat("\n===== RESUMO =====\n"); print(resumo)

#====Figuras====
paleta_div <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
paleta_seq <- hcl.colors(100, "YlOrRd", rev = TRUE)

mapa <- function(M, x, y, titulo, col = paleta_seq, simetrico = FALSE, zlim = NULL) {
    if (simetrico) { zm <- max(abs(M), na.rm = TRUE); zlim <- c(-zm, zm) }
    if (is.null(zlim)) zlim <- range(M, na.rm = TRUE)
    layout(matrix(1:2, 1), widths = c(6, 1))
    par(mar = c(4, 4, 3, 1))
    image(x, y, M, col = col, zlim = zlim, xlab = "Longitude", ylab = "Latitude", main = titulo, asp = 1)
    par(mar = c(4, 0.5, 3, 3.5))
    nz <- length(col)
    image(x = 0:1, y = seq(zlim[1], zlim[2], length.out = nz + 1),
          z = matrix(seq(zlim[1], zlim[2], length.out = nz), nrow = 1),
          col = col, zlim = zlim, axes = FALSE, xlab = "", ylab = "")
    axis(4); box()
    layout(1); par(mar = c(5, 4, 4, 2) + 0.1)
}

pdf(ARQ_PDF, width = 9.5, height = 7)

# --- Temporal: visão geral
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3, 1))
anos_ev <- ev$time / 365.25
plot(anos_ev, seq_len(n), type = "l", lwd = 2, xlab = "Anos desde o início", ylab = "Eventos acumulados",
     main = "N(t) observado vs. compensador")
lines(anos_ev, tempo_fs$tau, col = 2, lwd = 2, lty = 2)
legend("topleft", c("Observado N(t)", expression(hat(Lambda)(t))), col = c(1, 2), lty = c(1, 2), lwd = 2, bty = "n")
plot(tempo_fs$tau, seq_len(n), type = "l", lwd = 2, xlab = expression(tau[j]), ylab = "j",
     main = "Tempos transformados (esperado: diagonal)")
abline(0, 1, col = 2, lty = 2)
lines(tempo_fs$tau, tempo_fs$tau + 1.36 * sqrt(n), col = "grey50", lty = 3)
lines(tempo_fs$tau, tempo_fs$tau - 1.36 * sqrt(n), col = "grey50", lty = 3)
Uo <- sort(tempo_fs$U)
plot(ppoints(n), Uo, pch = 20, cex = 0.4, xlab = "Quantis U(0,1)", ylab = expression(U[j] ~ ordenado),
     main = sprintf("QQ-plot de U_j (KS: D=%.3f, p=%.3g)", tempo_fs$ks$statistic, tempo_fs$ks$p.value))
abline(0, 1, col = 2)
acf(tempo_fs$U, lag.max = 30, main = sprintf("ACF de U_j (Ljung-Box p=%.3g)", tempo_fs$lb$p.value))

# --- Temporal: janelas
par(mfrow = c(2, 1), mar = c(4.2, 4.2, 3, 1))
plot(anos_centro, O_t, type = "h", lwd = 4, col = "grey60", xlab = "Anos desde o início",
     ylab = "Eventos por janela", main = "Observado (barras) vs. esperado (linha)")
lines(anos_centro, E_t, col = 2, lwd = 2)
plot(anos_centro, R_t, type = "h", lwd = 4, col = ifelse(abs(R_t) > 2, 2, "grey60"),
     xlab = "Anos desde o início", ylab = "Resíduo de Pearson",
     main = "Resíduos de Pearson por janela (|R| > 2 em vermelho)")
abline(h = c(-2, 0, 2), lty = c(3, 1, 3))

# --- Comparação F_S real vs. F_S = 1
if (COMPARAR_SEM_FS) {
    par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))
    plot(tempo_fs$tau, tempo_fs$tau - seq_len(n), type = "l", lwd = 2, xlab = expression(tau[j]),
         ylab = expression(tau[j] - j), main = "Discrepância do compensador")
    lines(tempo_sem$tau, tempo_sem$tau - seq_len(n), col = 2, lwd = 2, lty = 2)
    abline(h = 0, lty = 3)
    legend("topleft", c("com F_S (exato)", "F_S = 1"), col = c(1, 2), lty = c(1, 2), lwd = 2, bty = "n")
    hist(FS, breaks = 30, col = "grey80", main = "Fração do kernel dentro da janela (F_S)",
         xlab = expression(F[S](i)))
}

# --- Espacial
par(mfrow = c(1, 1))
mapa(log10(O_ag + 1), x_ag, y_ag, sprintf("Observado: log10(contagem + 1), pixels %dx", AGREGAR))
mapa(log10(E_ag + 1), x_ag, y_ag, "Esperado: log10(contagem + 1)")
mapa(R_pix, x_ag, y_ag, sprintf("Resíduo de Pearson por pixel (E >= %g)", E_MIN),
     col = paleta_div, simetrico = TRUE)
mapa(res_suav, x_grid, y_grid, "Resíduo bruto suavizado (O - E)", col = paleta_div, simetrico = TRUE)
mapa(razao_log, x_grid, y_grid, "log[(O + 0.5)/(E + 0.5)] suavizado (E suavizado >= 0.5)",
     col = paleta_div, simetrico = TRUE)

par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3, 1))
hist(R_pix[!is.na(R_pix)], breaks = 30, freq = FALSE, col = "grey80", main = "Resíduos de Pearson espaciais",
     xlab = "R")
curve(dnorm(x), add = TRUE, col = 2, lwd = 2)
qqnorm(R_pix[!is.na(R_pix)], main = "QQ-plot normal (heurístico)", pch = 20, cex = 0.6); qqline(R_pix[!is.na(R_pix)], col = 2)
plot(x_ag, rowSums(O_ag), type = "h", lwd = 4, col = "grey60", xlab = "Longitude", ylab = "Eventos",
     main = "Perfil em longitude")
lines(x_ag, rowSums(E_ag), col = 2, lwd = 2)
plot(y_ag, colSums(O_ag), type = "h", lwd = 4, col = "grey60", xlab = "Latitude", ylab = "Eventos",
     main = "Perfil em latitude")
lines(y_ag, colSums(E_ag), col = 2, lwd = 2)

# --- Efeito de borda
par(mfrow = c(1, 2), mar = c(5, 4.2, 3, 1))
barplot(rho_faixa, names.arg = levels(faixa), col = "grey70", xlab = "Distância à borda da janela (graus)",
        ylab = "rho_j médio (prob. de background)", main = "rho_j vs. distância à borda", las = 2)
abline(h = mean(ev$rho_bg), col = 2, lty = 2)
par(mfrow = c(1, 1))
mapa(rho_ag, x_ag, y_ag, "rho_j médio por pixel (prob. de background)", zlim = c(0, 1))

dev.off()

#====Salvar====
res_peru <- list(
    resumo = resumo, FS = FS, tempo_fs = tempo_fs,
    tempo_sem = if (COMPARAR_SEM_FS) tempo_sem else NULL,
    tempo_bins = data.frame(centro_anos = anos_centro, O = O_t, E = E_t, R = R_t),
    espacial = list(O_pix = O_pix, E_pix = E_pix, O_ag = O_ag, E_ag = E_ag, R_pix = R_pix,
                    res_suav = res_suav, x_grid = x_grid, y_grid = y_grid, x_ag = x_ag, y_ag = y_ag),
    borda = list(d_borda = d_borda, rho_por_faixa = rho_faixa, n_por_faixa = n_faixa),
    config = list(M0 = M0, AGREGAR = AGREGAR, E_MIN = E_MIN, K_RAD = K_RAD, L_ANG = L_ANG,
                  SIGMA_SUAV = SIGMA_SUAV, SEMENTE = SEMENTE)
)
save(res_peru, file = ARQ_SAIDA, compress = TRUE)
cat("\nFiguras em", ARQ_PDF, "| resultados em", ARQ_SAIDA, "\n")
