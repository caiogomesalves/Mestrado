# =============================================================================
# Figuras de resultados do ETAS espaço-temporal (Peru)
# a partir de ajuste_peru.RData (objeto aj_peru).
#
# Saída: pasta figuras_peru/ com cada figura em PDF (vetorial) e PNG (300 dpi),
#        e tabela_parametros.csv.
#
# Figuras
#  01 Taxa de background contínua (eventos/ano/10^4 km^2), com costa e grandes eventos
#  02 Perfis da taxa de background (por latitude e por longitude)
#  03 Declustering estocástico: eventos de background vs. disparados
#  04 Fração de background por pixel
#  05 Série temporal anual: background vs. disparados
#  06 Kernels estimados: Omori-Utsu, produtividade, raio de disparo, densidade radial
#  07 Lei de Gutenberg-Richter e valor-b usado na razão de ramificação
#  08 Convergência do EM (trajetórias de theta e variação do raster)
#  09 Intensidade condicional total antes e depois do maior evento
#
# Observações
#  - Grade: o .RData não guarda x_grid/y_grid; reconstruída como no ajuste (precisa
#    de Inst_Peru.xlsx e do pacote ETAS) ou aproximada pela caixa dos eventos.
#  - Costa: usa o pacote 'maps' se instalado (opcional).
#  - Interpolação bilinear do raster é só visual; o raster é a densidade média por célula.
#  - Os erros-padrão da tabela vêm da Hessiana do passo M e SUBESTIMAM a incerteza
#    (ver seção 3.6): use-os como referência, não como intervalo final.
#  - Escrito sem execução em R no ambiente de redação; cada figura roda dentro de
#    tryCatch, então uma falha não derruba as demais.
# =============================================================================

suppressPackageStartupMessages(library(parallel))

#====Configurações====
ARQ_RDATA      <- "ajuste_peru.RData"
ARQ_XLSX       <- "Inst_Peru.xlsx"
DIR_FIG        <- "figuras_peru"
M0             <- 4.5
MAG_DESTAQUE   <- 7.0      # grandes eventos destacados nos mapas
FATOR_INTERP   <- 4        # refinamento visual do raster (bilinear)
DECADAS_ESCALA <- 3        # décadas mostradas na escala logarítmica de taxas
LIMIAR_RHO     <- 0.5      # corte para chamar um evento de background
N_MIN_PIXEL    <- 10       # eventos mínimos por pixel agregado no mapa de fração
AGREGAR        <- 4
NF             <- 200      # resolução da grade fina do mapa de intensidade total
DIAS_APOS      <- 1        # dias após o maior evento (figura 09)

#====Carregar====
load(ARQ_RDATA)
stopifnot(exists("aj_peru"))
th <- aj_peru$theta; ev <- aj_peru$eventos; raster <- aj_peru$raster
T_total <- aj_peru$T_total; n <- nrow(ev); Nx <- nrow(raster); Ny <- ncol(raster)
A_par <- unname(th["A"]); alpha_par <- unname(th["alpha"]); c_par <- unname(th["c"])
p_par <- unname(th["p"]); D_par <- unname(th["D"]); q_par <- unname(th["q"]); gamma_par <- unname(th["gamma"])
stopifnot(!is.unsorted(ev$time))

#====Grade e (opcional) data de início====
DATA_INICIO <- NULL
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
        DATA_INICIO <<- min(peru$date)
        cat("Grade reconstruída a partir de", ARQ_XLSX, "\n")
    } else {
        rl <- range(ev$x); rt <- range(ev$y)
        warning("Inst_Peru.xlsx indisponível: grade aproximada pela caixa dos eventos; ",
                "o eixo do tempo ficará em anos desde o início.")
    }
    list(x = seq(rl[1], rl[2], length.out = Nx), y = seq(rt[1], rt[2], length.out = Ny))
}
gr <- reconstruir_grade(); x_grid <- gr$x; y_grid <- gr$y
dx <- diff(x_grid)[1]; dy <- diff(y_grid)[1]
bx <- c(x_grid - dx / 2, x_grid[Nx] + dx / 2); by <- c(y_grid - dy / 2, y_grid[Ny] + dy / 2)
janela <- list(x0 = bx[1], x1 = bx[Nx + 1], y0 = by[1], y1 = by[Ny + 1])

#====Funções auxiliares====
kappa_func <- function(m, m0, A, alpha) A * exp(alpha * (m - m0))
g_dens <- function(t, c, p) (p - 1) * c^(p - 1) * (t + c)^(-p)
G_omori <- function(s, c, p) ifelse(s > 0, 1 - (c / (s + c))^(p - 1), 0)
S_de_m <- function(m) D_par * exp(gamma_par * (m - M0))

# Área (km^2) de 1 grau^2 em cada latitude; converte eventos/dia/grau^2 -> eventos/ano/10^4 km^2
area_km2_grau2 <- function(lat) 111.32 * cos(lat * pi / 180) * 110.57
conv_taxa <- function(M, y) M * 365.25 * 1e4 / matrix(area_km2_grau2(y), nrow(M), ncol(M), byrow = TRUE)

# Interpolação bilinear (tensorial) para visualização
bilinear <- function(M, x, y, xn, yn) {
    M1 <- apply(M, 2, function(col) approx(x, col, xn, rule = 2)$y)    # length(xn) x Ny
    t(apply(M1, 1, function(lin) approx(y, lin, yn, rule = 2)$y))      # length(xn) x length(yn)
}

desenhar_costa <- function() {
    if (requireNamespace("maps", quietly = TRUE)) {
        maps::map("world", xlim = c(janela$x0, janela$x1), ylim = c(janela$y0, janela$y1),
                  add = TRUE, col = "grey25", lwd = 0.8)
    }
}
grandes_eventos <- function(col = "black") {
    g <- ev$mag >= MAG_DESTAQUE
    if (any(g)) points(ev$x[g], ev$y[g], pch = 1, cex = 0.8 + (ev$mag[g] - MAG_DESTAQUE) * 1.2, col = col, lwd = 1.4)
}

paleta_taxa <- hcl.colors(100, "YlOrRd", rev = TRUE)
paleta_div  <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)

# Mapa com barra de cores; escala = "log10" (valores positivos) ou "linear"
mapa <- function(M, x, y, titulo, col = paleta_taxa, escala = c("linear", "log10"),
                 zlim = NULL, rotulo = "", extras = function() NULL) {
    escala <- match.arg(escala)
    Z <- if (escala == "log10") log10(pmax(M, .Machine$double.xmin)) else M
    if (is.null(zlim)) zlim <- range(Z[is.finite(Z)], na.rm = TRUE)
    Z <- pmin(pmax(Z, zlim[1]), zlim[2])
    layout(matrix(1:2, 1), widths = c(6, 1))
    par(mar = c(4, 4, 3, 0.5))
    image(x, y, Z, col = col, zlim = zlim, xlab = "Longitude", ylab = "Latitude", main = titulo, asp = 1)
    desenhar_costa(); extras(); box()
    par(mar = c(4, 0.5, 3, 4.5))
    nz <- length(col)
    image(x = 0:1, y = seq(zlim[1], zlim[2], length.out = nz + 1),
          z = matrix(seq(zlim[1], zlim[2], length.out = nz), nrow = 1),
          col = col, zlim = zlim, axes = FALSE, xlab = "", ylab = "")
    if (escala == "log10") {
        tk <- ceiling(zlim[1]):floor(zlim[2])
        axis(4, at = tk, labels = parse(text = paste0("10^", tk)), las = 1)
    } else axis(4, las = 1)
    mtext(rotulo, side = 4, line = 3.2, cex = 0.8); box()
    layout(1); par(mar = c(5, 4, 4, 2) + 0.1)
}

salvar_fig <- function(nome, fun, w = 7, h = 6) {
    dir.create(DIR_FIG, showWarnings = FALSE)
    r <- tryCatch({
        pdf(file.path(DIR_FIG, paste0(nome, ".pdf")), width = w, height = h); fun(); dev.off()
        png(file.path(DIR_FIG, paste0(nome, ".png")), width = w, height = h, units = "in", res = 300); fun(); dev.off()
        TRUE
    }, error = function(e) { while (!is.null(dev.list())) dev.off(); message("Falha em ", nome, ": ", conditionMessage(e)); FALSE })
    if (isTRUE(r)) cat("ok:", nome, "\n")
}

#====Quantidades derivadas====
xf  <- seq(janela$x0, janela$x1, length.out = Nx * FATOR_INTERP)
yf  <- seq(janela$y0, janela$y1, length.out = Ny * FATOR_INTERP)
taxa_nativa <- conv_taxa(raster, y_grid)                          # eventos/ano/10^4 km^2
taxa_fina   <- conv_taxa(bilinear(raster, x_grid, y_grid, xf, yf), yf)

# valor-b (Aki-Utsu) e razão de ramificação
dm_bin <- min(diff(sort(unique(round(ev$mag, 2)))))
b_hat    <- log10(exp(1)) / (mean(ev$mag) - (M0 - dm_bin / 2))
beta_hat <- b_hat * log(10)
n_ramif  <- if (beta_hat > alpha_par) A_par * beta_hat / (beta_hat - alpha_par) else Inf

#====Tabela de parâmetros====
lower <- c(A = 1e-6, alpha = 1e-6, c = 1e-6, p = 1.001, D = 1e-6, q = 1.001, gamma = 1e-6)
upper <- c(A = 15, alpha = 10, c = 15, p = 15, D = 15, q = 15, gamma = 10)
H <- aj_peru$hessiana
covH <- tryCatch(solve(H), error = function(e) NULL)
se <- if (!is.null(covH) && all(diag(covH) > 0)) sqrt(diag(covH)) else rep(NA_real_, length(th))
names(se) <- names(th)
na_cota <- (th - lower[names(th)]) < 1e-3 * abs(lower[names(th)] + 1e-12) + 1e-6 |
           (upper[names(th)] - th) < 1e-3 * abs(upper[names(th)])
tab_par <- data.frame(parametro = names(th), estimativa = as.numeric(th), ep_hessiana = as.numeric(se),
                      ic95_inf = as.numeric(th) - 1.96 * as.numeric(se),
                      ic95_sup = as.numeric(th) + 1.96 * as.numeric(se),
                      sobre_cota = as.logical(na_cota))
write.csv(tab_par, "tabela_parametros.csv", row.names = FALSE)
cat("\n===== Parâmetros (ep da Hessiana do passo M: subestimam a incerteza) =====\n"); print(tab_par, digits = 4)

se_n <- NA_real_
if (!is.null(covH) && is.finite(n_ramif)) {
    g <- c(beta_hat / (beta_hat - alpha_par), A_par * beta_hat / (beta_hat - alpha_par)^2)
    se_n <- sqrt(drop(t(g) %*% covH[c("A", "alpha"), c("A", "alpha")] %*% g))
}
cat(sprintf("\nvalor-b (Aki-Utsu, bin %.2f) = %.3f | razão de ramificação n = %.3f (ep delta ~ %.3f; beta fixo)\n",
            dm_bin, b_hat, n_ramif, se_n))
cat(sprintf("Fração de background (média de rho_j) = %.3f | 1 - n = %.3f (iguais no regime estacionário)\n",
            mean(ev$rho_bg), 1 - n_ramif))
cat(sprintf("Convergiu: %s | iterações do EM: %d\n", aj_peru$convergiu, aj_peru$n_iter_em))

#====01 Taxa de background contínua====
salvar_fig("01_background_taxa", function() {
    mapa(taxa_fina, xf, yf, "Taxa de background estimada",
         escala = "log10", rotulo = "eventos / ano / 10^4 km^2",
         zlim = { zm <- ceiling(log10(max(taxa_fina))); c(zm - DECADAS_ESCALA, zm) },
         extras = function() grandes_eventos("black"))
}, w = 7.5, h = 6.5)

#====02 Perfis da taxa de background====
salvar_fig("02_background_perfis", function() {
    par(mfrow = c(1, 2), mar = c(4.2, 4.5, 3, 1))
    por_lat <- rowSums(t(raster)) * dx * 365.25              # eventos/ano por grau de latitude
    por_lon <- rowSums(raster) * dy * 365.25                 # eventos/ano por grau de longitude
    plot(por_lat, y_grid, type = "l", lwd = 2, xlab = "Eventos / ano / grau de latitude", ylab = "Latitude",
         main = "Perfil por latitude")
    polygon(c(0, por_lat, 0), c(y_grid[1], y_grid, y_grid[Ny]), col = adjustcolor("#D95F0E", 0.35), border = NA)
    plot(x_grid, por_lon, type = "l", lwd = 2, xlab = "Longitude", ylab = "Eventos / ano / grau de longitude",
         main = "Perfil por longitude")
    polygon(c(x_grid[1], x_grid, x_grid[Nx]), c(0, por_lon, 0), col = adjustcolor("#D95F0E", 0.35), border = NA)
}, w = 9, h = 4.5)

#====03 Declustering: background vs. disparados====
salvar_fig("03_declustering_mapa", function() {
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    bg <- ev$rho_bg >= LIMIAR_RHO
    cex_m <- 0.25 + (ev$mag - M0) * 0.35
    for (k in 1:2) {
        sel <- if (k == 1) bg else !bg
        plot(NA, xlim = c(janela$x0, janela$x1), ylim = c(janela$y0, janela$y1), asp = 1,
             xlab = "Longitude", ylab = "Latitude",
             main = sprintf("%s (n = %d; %.1f%%)", if (k == 1) "Background" else "Disparados", sum(sel), 100 * mean(sel)))
        desenhar_costa()
        points(ev$x[sel], ev$y[sel], pch = 16, cex = cex_m[sel],
               col = adjustcolor(if (k == 1) "#2166AC" else "#B2182B", 0.35))
        grandes_eventos("black"); box()
    }
}, w = 11, h = 6)

#====04 Fração de background por pixel====
salvar_fig("04_fracao_background_mapa", function() {
    ix <- findInterval(ev$x, bx, all.inside = TRUE); iy <- findInterval(ev$y, by, all.inside = TRUE)
    agr <- function(M, f) {
        nx <- floor(nrow(M) / f) * f; ny <- floor(ncol(M) / f) * f; M <- M[1:nx, 1:ny, drop = FALSE]
        out <- matrix(0, nx / f, ny / f)
        for (a in 1:f) for (b in 1:f) out <- out + M[seq(a, nx, by = f), seq(b, ny, by = f)]
        out
    }
    cont <- agr(matrix(tabulate(ix + (iy - 1) * Nx, nbins = Nx * Ny), Nx, Ny), AGREGAR)
    soma <- agr(matrix(as.numeric(tapply(ev$rho_bg, factor(ix + (iy - 1) * Nx, levels = 1:(Nx * Ny)), sum)),
                       Nx, Ny), AGREGAR)
    soma[is.na(soma)] <- 0
    fr <- ifelse(cont >= N_MIN_PIXEL, soma / cont, NA)
    xa <- colMeans(matrix(x_grid[1:(nrow(fr) * AGREGAR)], nrow = AGREGAR))
    ya <- colMeans(matrix(y_grid[1:(ncol(fr) * AGREGAR)], nrow = AGREGAR))
    mapa(fr, xa, ya, sprintf("Fração de background por pixel (n >= %d eventos)", N_MIN_PIXEL),
         col = hcl.colors(100, "Blues 3", rev = TRUE), zlim = c(0, 1), rotulo = "fração de background")
}, w = 7.5, h = 6.5)

#====05 Série temporal anual====
salvar_fig("05_serie_temporal", function() {
    if (!is.null(DATA_INICIO)) {
        ano <- as.integer(format(DATA_INICIO + ev$time, "%Y")); rot_x <- "Ano"
    } else {
        ano <- floor(ev$time / 365.25); rot_x <- "Anos desde o início do catálogo"
    }
    anos <- sort(unique(ano))
    bg <- tapply(ev$rho_bg, factor(ano, levels = anos), sum)
    tr <- tapply(1 - ev$rho_bg, factor(ano, levels = anos), sum)
    par(mar = c(5, 4.5, 3, 1))
    bp <- barplot(rbind(bg, tr), col = c("#2166AC", "#B2182B"), border = NA, space = 0.1,
                  names.arg = ifelse(anos %% 5 == 0, anos, ""), las = 2, xlab = rot_x,
                  ylab = "Eventos por ano (esperados pelos pesos)", main = "Eventos anuais: background e disparados")
    legend("topleft", c("Background (soma de rho)", "Disparados (soma de 1 - rho)"),
           fill = c("#2166AC", "#B2182B"), bty = "n")
}, w = 9, h = 4.8)

#====06 Kernels estimados====
salvar_fig("06_kernels", function() {
    par(mfrow = c(2, 2), mar = c(4.2, 4.5, 3, 1))
    # (a) Omori-Utsu
    tt <- 10^seq(-3, 4, length.out = 300)
    plot(tt, g_dens(tt, c_par, p_par), type = "l", log = "xy", lwd = 2, xlab = "Tempo após o evento (dias)",
         ylab = "g(t)", main = sprintf("Omori-Utsu (c = %.3g, p = %.3f)", c_par, p_par))
    abline(v = c_par, lty = 3, col = "grey50")
    legend("bottomleft", sprintf("fração em 1 dia: %.2f | 30 dias: %.2f | 1 ano: %.2f",
           G_omori(1, c_par, p_par), G_omori(30, c_par, p_par), G_omori(365, c_par, p_par)),
           bty = "n", cex = 0.75)
    # (b) produtividade
    mm <- seq(M0, 8.5, length.out = 100)
    plot(mm, kappa_func(mm, M0, A_par, alpha_par), type = "l", log = "y", lwd = 2, xlab = "Magnitude",
         ylab = "Filhos diretos esperados", main = sprintf("Produtividade (A = %.3g, alpha = %.2f)", A_par, alpha_par))
    # (c) raio de disparo
    rq <- function(m, prob) sqrt(S_de_m(m) * ((1 - prob)^(-1 / (q_par - 1)) - 1)) * 111
    plot(mm, rq(mm, 0.9), type = "l", log = "y", lwd = 2, lty = 2, xlab = "Magnitude", ylab = "Raio (km)",
         main = sprintf("Raio de disparo (q = %.2f)", q_par), ylim = range(c(rq(mm, 0.5), rq(mm, 0.9))))
    lines(mm, rq(mm, 0.5), lwd = 2)
    legend("topleft", c("50% dos filhos", "90% dos filhos"), lty = c(1, 2), lwd = 2, bty = "n")
    # (d) densidade radial
    rr <- 10^seq(-1, 3, length.out = 300)
    plot(NA, xlim = range(rr), ylim = c(1e-6, 1), log = "xy", xlab = "Distância (km)", ylab = "Densidade radial (por km)",
         main = "Densidade radial do kernel de Pareto")
    cols <- c("#2166AC", "#1B9E77", "#B2182B"); ms <- c(5, 6, 7)
    for (k in 1:3) {
        S <- S_de_m(ms[k]); r_deg <- rr / 111
        dens <- 2 * r_deg * (q_par - 1) / S * (1 + r_deg^2 / S)^(-q_par) / 111
        lines(rr, dens, col = cols[k], lwd = 2)
    }
    legend("bottomleft", paste("M =", ms), col = cols, lwd = 2, bty = "n")
}, w = 9, h = 7)

#====07 Gutenberg-Richter====
salvar_fig("07_gutenberg_richter", function() {
    ms <- seq(M0, max(ev$mag), by = dm_bin)
    N_acum <- vapply(ms, function(m) sum(ev$mag >= m - 1e-9), numeric(1))
    par(mar = c(4.5, 4.5, 3, 1))
    plot(ms, N_acum, log = "y", pch = 16, cex = 0.7, xlab = "Magnitude (M)", ylab = "N(>= M)",
         main = sprintf("Gutenberg-Richter: b = %.2f (Aki-Utsu)", b_hat))
    lines(ms, n * 10^(-b_hat * (ms - M0)), col = 2, lwd = 2)
    legend("topright", c("Observado", "Ajuste"), pch = c(16, NA), lty = c(NA, 1), col = c(1, 2), bty = "n")
}, w = 6.5, h = 5)

#====08 Convergência do EM====
salvar_fig("08_convergencia_EM", function() {
    if (is.null(aj_peru$thetas_hist)) stop("thetas_hist ausente (ajuste com guardar_historico = TRUE).")
    TH <- do.call(rbind, aj_peru$thetas_hist)
    par(mfrow = c(3, 3), mar = c(4, 4.2, 2.5, 1))
    for (k in seq_len(ncol(TH))) {
        plot(0:(nrow(TH) - 1), TH[, k], type = "b", pch = 16, cex = 0.6, xlab = "Iteração", ylab = colnames(TH)[k],
             main = colnames(TH)[k])
    }
    if (!is.null(aj_peru$rasters_hist)) {
        RH <- aj_peru$rasters_hist
        dif <- vapply(2:length(RH), function(i) max(abs(RH[[i]] - RH[[i - 1]])) / max(abs(RH[[i - 1]])), numeric(1))
        plot(seq_along(dif), dif, type = "b", pch = 16, cex = 0.6, log = "y", xlab = "Iteração",
             ylab = "max |Delta raster| relativo", main = "Variação do raster")
    }
}, w = 10, h = 8)

#====09 Intensidade condicional total antes e depois do maior evento====
salvar_fig("09_intensidade_total", function() {
    i_ref <- which.max(ev$mag)
    xi <- seq(janela$x0, janela$x1, length.out = NF); yi <- seq(janela$y0, janela$y1, length.out = NF)
    mu_f <- bilinear(raster, x_grid, y_grid, xi, yi)
    kap_i <- kappa_func(ev$mag, M0, A_par, alpha_par); S_all <- S_de_m(ev$mag)
    lambda_em <- function(t) {
        idx <- which(ev$time < t)
        peso <- kap_i[idx] * g_dens(t - ev$time[idx], c_par, p_par)
        proxy <- peso / S_all[idx]                         # ~ intensidade de pico de cada contribuição
        manter <- idx[proxy > 1e-6 * max(proxy)]
        L <- mu_f
        for (k in manter) {
            d2 <- outer((xi - ev$x[k])^2, (yi - ev$y[k])^2, "+")
            L <- L + kap_i[k] * g_dens(t - ev$time[k], c_par, p_par) *
                (q_par - 1) / (pi * S_all[k]) * (1 + d2 / S_all[k])^(-q_par)
        }
        conv_taxa(L, yi)
    }
    t_ref <- ev$time[i_ref]
    Lant <- lambda_em(t_ref - 1e-6); Lpos <- lambda_em(t_ref + DIAS_APOS)
    zm <- ceiling(log10(max(Lpos))); zlim <- c(zm - 2 * DECADAS_ESCALA, zm)
    par(mfrow = c(1, 1))
    mapa(Lant, xi, yi, sprintf("Intensidade logo antes do maior evento (M %.1f)", ev$mag[i_ref]), escala = "log10",
         zlim = zlim, rotulo = "eventos / ano / 10^4 km^2",
         extras = function() points(ev$x[i_ref], ev$y[i_ref], pch = 8, col = "black", cex = 1.5, lwd = 1.5))
    mapa(Lpos, xi, yi, sprintf("Intensidade %g dia(s) depois", DIAS_APOS), escala = "log10",
         zlim = zlim, rotulo = "eventos / ano / 10^4 km^2",
         extras = function() points(ev$x[i_ref], ev$y[i_ref], pch = 8, col = "black", cex = 1.5, lwd = 1.5))
}, w = 7.5, h = 6.5)

cat("\nFiguras em", DIR_FIG, "| tabela em tabela_parametros.csv\n")
