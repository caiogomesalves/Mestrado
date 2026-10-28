#----Pacotes----

library(tidyverse)
library(ETAS)
library(spatstat)
library(readxl)
library(parallel)
library(httr)

#----Ajuste----

region_iran <- list(lat = c(26, 25, 29, 38, 35),
                    long = c(52, 59, 58, 45, 43))

# Catálogo:
iran.cat <- catalog(iran.quakes,
                    study.start = "1991/01/01", study.end = "2011/01/01",
                    region.poly = region_iran,
                    mag.threshold = 4.5)

plot(iran.cat)

# Chute inicial parâmetros:
theta_iran <- c("A" = 0.2, "alpha" = 2.7, "c" = 0.02,
                "p" = 1.2, "D" = 0.02, "q" = 2.3, "gamma" = 0.03)

thetas_iran <- list()
thetas_iran[[1]] <- theta_iran

# Número total de eventos:
n_events_iran <- nrow(iran.cat$revents[iran.cat$revents[, "flag"] %in% c(0, 1), ])

# Definir as sequências do grid explicitamente
x_grid_iran <- seq(min(iran.cat$longlat.coord$long), max(iran.cat$longlat.coord$long), length.out = 64)
y_grid_iran <- seq(min(iran.cat$longlat.coord$lat), max(iran.cat$longlat.coord$lat), length.out = 64)

# Raster inicial:
r_iran <- matrix(1, nrow = 128, ncol = 128)
rel_clust_iran <- matrix(0, nrow = 128, ncol = 128)

# Normaliza o raster para a integral ser 1:
# r_iran <- r_iran/(sum(r_iran))

# Na lista de rasters, guardaremos as matrizes:
rasters_iran <- list()
rast_rel_clust_iran <- list()
rasters_iran[[1]] <- r_iran

# Data.frame com os eventos para funções:
df_temp_iran <- data.frame(x = iran.cat$longlat.coord$long[iran.cat$revents[, "flag"] %in% c(0, 1)],
                           y = iran.cat$longlat.coord$lat[iran.cat$revents[, "flag"] %in% c(0, 1)],
                           mag = iran.cat$revents[iran.cat$revents[, "flag"] %in% c(0, 1), "mm"] + 4.5,
                           time = iran.cat$revents[iran.cat$revents[, "flag"] %in% c(0, 1), "tt"])

df_temp_iran$time <- df_temp_iran$time - min(df_temp_iran$time)

# Distâncias para kernel adaptativo:
d_j_iran <- d_j(iran.cat$longlat.coord[iran.cat$revents[, "flag"] %in% c(0, 1), c("long", "lat")])

num_cores <- 8

t1_iran <- Sys.time()

# Paralelização do cálculo de probabilidades:
probabilidades_iran <- mclapply(1:n_events_iran, function(j) {
    # Primeiro evento tem probabilidade 0 de ser descendente:
    if (j == 1) return(numeric(0))
    # Vetor de índices dos antecessores:
    i_indices <- 1:(j - 1)
    # Denominador comum para cada evento j:
    denom <- lambda_total(
        df_temp_iran$time[j], df_temp_iran$x[j], df_temp_iran$y[j], df_temp_iran$mag[j],
        4.5, df_temp_iran, theta_iran, r_iran, x_grid_iran, y_grid_iran
    )
    # Computação do numerador para todos os antecessores:
    numeradores <- kappa_func(df_temp_iran$mag[i_indices], 4.5, theta_iran["A"], theta_iran["alpha"]) *
        g_func(df_temp_iran$time[j] - df_temp_iran$time[i_indices], theta_iran["c"], theta_iran["p"]) *
        f_spatial(
        (df_temp_iran$x[j] - df_temp_iran$x[i_indices])^2 + (df_temp_iran$y[j] - df_temp_iran$y[i_indices])^2,
        df_temp_iran$mag[i_indices], 4.5, theta_iran["D"], theta_iran["q"], theta_iran["gamma"]
        )
    # Retorna o vetor p_{i,j} para o evento j:
    return(numeradores / denom)
}, mc.cores = num_cores)

# Probabilidades totais:
prob_total_iran <- probabilidades_iran %>%
    lapply(sum) %>%
    do.call(c, .)

# Estima o background rate:
r_iran <- estimate_background_kernel(x_grid_iran,
                                     y_grid_iran,
                                     data.frame(x = iran.cat$longlat.coord$long[iran.cat$revents[, "flag"] %in% c(0, 1)],
                                                y = iran.cat$longlat.coord$lat[iran.cat$revents[, "flag"] %in% c(0, 1)],
                                                rho = prob_total_iran),
                                     d_j_iran,
                                     T_total = max(df_temp_iran$time))

rel_clust_iran <- estimate_relative_clustering(x_grid_iran,
                                               y_grid_iran,
                                               data.frame(x = iran.cat$longlat.coord$long[iran.cat$revents[, "flag"] %in% c(0, 1)],
                                                          y = iran.cat$longlat.coord$lat[iran.cat$revents[, "flag"] %in% c(0, 1)],
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
                         lower = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-6),
                         upper = c(15, 10, 15, 15, 15, 15, 10),
                         control = list(maxit = 1000),
                         hessian = T)

theta_iran <- otimizacao_iran$par
thetas_iran[[2]] <- theta_iran

rast_rel_clust_iran[[1]] <- rel_clust_iran

rasters_iran[[2]] <- r_iran

i <- 2

# Considerando tolerância de 10e-6:
while ((any(abs(rasters_iran[[i]] - rasters_iran[[i - 1]]) > 1e-3)) | (any(abs(thetas_iran[[i]] - thetas_iran[[i - 1]]) > 1e-3))) {
    # Calcula pares de probabilidades:
    probabilidades_iran <- mclapply(1:n_events_iran, function(j) {
        if (j == 1) return(numeric(0))
        i_indices <- 1:(j - 1)
        denom <- lambda_total(
            df_temp_iran$time[j], df_temp_iran$x[j], df_temp_iran$y[j], df_temp_iran$mag[j],
            4.5, df_temp_iran, theta_iran, r_iran, x_grid_iran, y_grid_iran
        )
        numeradores <- kappa_func(df_temp_iran$mag[i_indices], 4.5, theta_iran["A"], theta_iran["alpha"]) *
            g_func(df_temp_iran$time[j] - df_temp_iran$time[i_indices], theta_iran["c"], theta_iran["p"]) *
            f_spatial(
            (df_temp_iran$x[j] - df_temp_iran$x[i_indices])^2 + (df_temp_iran$y[j] - df_temp_iran$y[i_indices])^2,
            df_temp_iran$mag[i_indices], 4.5, theta_iran["D"], theta_iran["q"], theta_iran["gamma"]
            )
        return(numeradores / denom)
    }, mc.cores = num_cores)
    # Calcula probabilidades totais:
    prob_total_iran <- probabilidades_iran %>%
        lapply(sum) %>%
        do.call(c, .)
    # Estima parâmetros theta:
    otimizacao_iran <- optim(theta_iran,
                             log_lik_etas,
                             gr = gradiente_etas_completo,
                             df_eventos = df_temp_iran,
                             matriz_p_ij = probabilidades_iran,
                             prob_total = prob_total_iran,
                             raster_mu = r_iran, x_grid = x_grid_iran, y_grid = y_grid_iran,
                             method = "L-BFGS-B",
                             lower = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-6),
                             upper = c(15, 10, 15, 15, 15, 15, 10),
                             control = list(maxit = 1000),
                             hessian = T)
    theta_iran <- otimizacao_iran$par
    thetas_iran[[i + 1]] <- theta_iran
    # Estima background rate:
    r_iran <- estimate_background_kernel(x_grid_iran,
                                         y_grid_iran,
                                         data.frame(x = iran.cat$longlat.coord$long[iran.cat$revents[, "flag"] %in% c(0, 1)],
                                                    y = iran.cat$longlat.coord$lat[iran.cat$revents[, "flag"] %in% c(0, 1)],
                                                    rho = prob_total_iran),
                                         d_j_iran,
                                         T_total = max(df_temp_iran$time))
    rasters_iran[[i + 1]] <- r_iran
    # Estima clustering relativo:
    rel_clust_iran <- estimate_relative_clustering(x_grid_iran,
                                                   y_grid_iran,
                                                   data.frame(x = iran.cat$longlat.coord$long[iran.cat$revents[, "flag"] %in% c(0, 1)],
                                                              y = iran.cat$longlat.coord$lat[iran.cat$revents[, "flag"] %in% c(0, 1)],
                                                              rho = prob_total_iran),
                                                   d_j_iran)$C
    rast_rel_clust_iran[[i]] <- rel_clust_iran
    # Mostra progresso do loop:
    print(i + 1)
    # Atualiza índice:
    i <- i + 1
}

#stopCluster(cl)

t2_iran <- Sys.time()

#====Ajuste com função única====
t1_iran <- Sys.time()

aj_iran <- ajustar_etas_proprio(catalogo = iran.cat, theta0 = theta_iran, M0 = 4.5,
                                x_grid = x_grid_iran, y_grid = y_grid_iran, num_cores = 8)

t2_iran <- Sys.time()

aj_iran_etas <- etas(iran.cat, param0 = c(0.5, 0.2, 0.05, 2.7, 1.2, 0.02, 2.3, 0.03),
                     nthreads = 8)

t3_iran <- Sys.time()

#====Análise de Resíduos====

# Taxa global de background (Eventos/dia em toda a região):
mu_global_iran <- sum(pmax(0, 1 - prob_total_iran)) / max(df_temp_iran$time)

# Tau_j para os tempos dos eventos:
tau_j_iran <- calcular_tau(df_temp_iran$time, df_temp_iran, thetas_iran[[length(thetas_iran)]], mu_global_iran, M0 = 4.5)

# Calcular U_j:
delta_tau_iran <- diff(c(0, tau_j_iran))
U_j_iran <- 1 - exp(-delta_tau_iran)

#----Gráficos----
par(mfrow = c(1, 2))

# Gráfico A: tau_j vs j
plot(1:length(tau_j_iran), tau_j_iran, type = "l", col = "blue", lwd = 2,
     xlab = "Número do Evento (j)", ylab = expression(tau[j]),
     main = "Tempos Transformados")
abline(a = 0, b = 1, col = "red", lty = 2)

# Gráfico B: Q-Q Plot de U_j
plot(sort(U_j_iran), qunif(ppoints(length(U_j_iran))), pch = 20, col = "darkgreen",
     xlab = "Empírico", ylab = "Teórico (Uniforme 0-1)",
     main = "Q-Q Plot de U_j")
abline(a = 0, b = 1, col = "red", lty = 2)

ks.test(U_j_iran, punif)

# Número de caixas (bins) de tempo
n_temp_iran <- 100
janelas_tempo_iran <- seq(0, max(df_temp_iran$time), length.out = n_temp_iran + 1)

# Contagem de eventos em cada janela:
obs_temp_iran <- hist(df_temp_iran$time, breaks = janelas_tempo_iran, plot = FALSE)$counts

# Contagens esperadas pela diferença do tau nos limites da janela:
tau_limites_iran <- calcular_tau(janelas_tempo_iran, df_temp_iran, thetas_iran[[length(thetas_iran)]], mu_global_iran, M0 = 4.5)
esp_temp_iran <- diff(tau_limites_iran)

# Resíduos de Pearson Temporais:
res_temp_pearson_iran <- (obs_temp_iran - esp_temp_iran) / sqrt(esp_temp_iran)

# Resíduos Temporais:
par(mfrow = c(1, 1))
plot(janelas_tempo_iran[-1], res_temp_pearson_iran, type = "h", lwd = 2, col = "purple",
     xlab = "Tempo (dias)", ylab = "Resíduos de Pearson",
     main = "Resíduos Temporais Suavizados")
abline(h = 0, col = "black")
abline(h = c(-2, 2), col = "red", lty = 2)

#====Espaciais====

# Definir os pontos observados:
pontos_iran <- ppp(df_temp_iran$x, df_temp_iran$y,
                   window = owin(x_grid_iran[c(1, length(x_grid_iran))],
                                 y_grid_iran[c(1, length(y_grid_iran))]),
              marks = pmax(0, 1 - prob_total_iran))

# Criar a densidade espacial esperada pelo Background:
cell_area_iran <- diff(x_grid_iran)[1] * diff(y_grid_iran)[1]
matriz_esperada_iran <- (rasters_iran[[length(rasters_iran)]] / (sum(rasters_iran[[length(rasters_iran)]]) * cell_area_iran)) * (mu_global_iran * max(df_temp_iran$time))
# matriz_esperada_iran <- (rasters_iran[[length(rasters_iran)]] / sum(rasters_iran[[length(rasters_iran)]])) * (mu_global_iran * max(df_temp_iran$time))

# Extrair os 7 parâmetros do modelo convergido
A_iran <- thetas_iran[[length(thetas_iran)]]["A"]
alpha_iran <- thetas_iran[[length(thetas_iran)]]["alpha"]
D_iran <- thetas_iran[[length(thetas_iran)]]["D"]
q_par_iran <- thetas_iran[[length(thetas_iran)]]["q"]
gamma_par_iran <- thetas_iran[[length(thetas_iran)]]["gamma"]

M0 <- 4.5 # Ajuste para o limiar de magnitude correto do seu catálogo atual

# Sobrepor a nuvem espacial de Pareto (Lei de Potência) de TODOS os eventos
for(i in 1:nrow(df_temp_iran)) {
    xi_iran <- df_temp_iran$x[i]
    yi_iran <- df_temp_iran$y[i]
    mi_iran <- df_temp_iran$mag[i]
    # Produtividade esperada do evento (mantém-se com alpha)
    k_i_iran <- A_iran * exp(alpha_iran * (mi_iran - M0))
    # Fator de escala espacial S_i (agora independente de alpha, usa gamma)
    S_i_iran <- D_iran * exp(gamma_par_iran * (mi_iran - M0))
    # Matriz de distâncias ao quadrado
    dx2_iran <- (x_grid_iran - xi_iran)^2
    dy2_iran <- (y_grid_iran - yi_iran)^2
    dist2_matrix_iran <- outer(dx2_iran, dy2_iran, "+")
    # Kernel de Pareto Bivariado Isotrópico (Fórmula do pacote ETAS)
    # f(x,y) = ((q-1) / (pi * S_i)) * (1 + dist2 / S_i)^(-q)
    term1_iran <- (q_par_iran - 1) / (pi * S_i_iran)
    term2_iran <- (1 + (dist2_matrix_iran / S_i_iran))^(-q_par_iran)
    kernel_esp_iran <- k_i_iran * term1_iran * term2_iran
    matriz_esperada_iran <- matriz_esperada_iran + kernel_esp_iran
}

densidade_esp_iran <- im(t(matriz_esperada_iran), xcol = x_grid_iran, yrow = y_grid_iran)

# Criar a densidade observada, forçando o alinhamento com xy:
densidade_obs_iran <- density(pontos_iran, sigma = bw.diggle(pontos_iran), xy = densidade_esp_iran,
                              weights = pontos_iran$marks)

# Calcular os resíduos de Pearson:
epsilon <- 1e-10
res_spat_pearson_iran <- eval.im((densidade_obs_iran - densidade_esp_iran) / sqrt(densidade_esp_iran + epsilon))

# Visualização:

min_val <- min(res_spat_pearson_iran)
max_val <- max(res_spat_pearson_iran)

max_abs <- max(abs(c(min_val, max_val)))
sym_breaks <- seq(-max_abs, max_abs, length.out = 257)

# 2. Build a dedicated colourmap object mapping the breaks to the Oslo palette
col_map <- colourmap(
    col = hcl.colors(256, palette = "Viridis", rev = F),
  breaks = sym_breaks
)

plot(res_spat_pearson_iran, main = "Resíduos Espaciais Totais de Pearson (Modelo Pareto)", col = col_map)
plot(res_spat_pearson_iran, main = "Resíduos Espaciais Totais de Pearson (Modelo Pareto)", col = topo.colors(256))

#====Tolerância====

diferencas_rasters_iran <- list()

for (i in 2:length(rasters_iran)) {
    diferencas_rasters_iran[[i]] <- max(abs(rasters_iran[[i]] - rasters_iran[[i - 1]]))
}

diferencas_thetas_iran <- list()

for (i in 2:length(thetas_iran)) {
    diferencas_thetas_iran[[i]] <- max(abs(thetas_iran[[i]] - thetas_iran[[i - 1]]))
}

do.call(c, diferencas_rasters_iran) %>%
    cbind(do.call(c, diferencas_thetas_iran)) %>%
    matplot(type = "b", pch = c(1, 1))
legend("topright", legend = c("Rasteres", "Parâmetros"),
       col = c("black", "red"), pch = c(1, 1))
abline(h = 1e-3)

ts(do.call(rbind, thetas_iran)) %>%
    lattice::xyplot(xlab = "Iteração", pch = 16, type = "b",
                    scales = list(x = list(tick.number = length(thetas)),
                                  y = list(tick.number = 3)),
                    main = "Estimativa dos parâmetros do modelo")

#====Erros-Padrão====

otimizacao_iran$hessian %>% solve() %>% diag() %>% sqrt()

#====Ajuste pacote ETAS====

param01_iran <- c(0.5, 0.2, 0.05, 2.7, 1.2, 0.02, 2.3, 0.03)

iran.fit <- etas(iran.cat, param0 = param01_iran)

plot(iran.fit)

im(t(rasters_iran[[15]]), xrange = long.range, yrange = lat.range) %>% plot()

#====Diagnóstico visual do background rate====

grid_iran <- expand.grid(
    x = x_grid_iran,
    y = y_grid_iran
)

df_iran <- extrair_eventos_catalog(iran.cat, 4.5)

par(mfrow = c(2, 1), mar = c(4, 4, 2, 5))

plot(df_iran$x, df_iran$y, xlab = "long",
     ylab = "lat", col = "blue", cex = 2 * (df_iran$mag - 4.4)/max(df_iran$mag - 4.4),
     asp = TRUE, axes = FALSE,
     xlim = range(df_iran$x),
     ylim = range(df_iran$y))
maps::map("world", add = TRUE, col = "grey50")
axis(1)
axis(2)
box()
fields::quilt.plot(grid_iran$x, grid_iran$y, aj_iran$raster,
           nx = 64, ny = 64, asp = T,
           xlim = range(df_iran$x),
           ylim = range(df_iran$y),
           add.legend = F)
maps::map("world", add = TRUE, col = "grey50")
axis(1)
axis(2)
box()
fields::image.plot(zlim = range(aj_iran$raster, na.rm = TRUE),
                   legend.only = TRUE)

#----Perfis de verossimilhança----

grade_A <- seq(1e-6, 1, length.out = 50)
grade_alpha <- seq(1, 5, length.out = 50)
grade_c <- seq(1e-6, 1, length.out = 50)
grade_p <- seq(1 + 1e-6, 2, length.out = 50)
grade_D <- seq(1e-6, 1, length.out = 50)
grade_q <- seq(1 + 1e-6, 5, length.out = 50)
grade_gamma <- seq(1e-6, 1, length.out = 50)

perfil_A <- perfilar_rapido(
    ajuste = aj_iran, df_eventos = df_iran,
    M0 = 4.5,
    x_grid = x_grid_iran, y_grid = y_grid_iran,
    parametro = "A", grade_valores = grade_A,
    num_cores = 8
)

perfil_alpha <- perfilar_rapido(
    ajuste = aj_iran, df_eventos = df_iran,
    M0 = 4.5,
    x_grid = x_grid_iran, y_grid = y_grid_iran,
    parametro = "alpha", grade_valores = grade_alpha,
    num_cores = 8
)

perfil_c <- perfilar_rapido(
    ajuste = aj_iran, df_eventos = df_iran,
    M0 = 4.5,
    x_grid = x_grid_iran, y_grid = y_grid_iran,
    parametro = "c", grade_valores = grade_c,
    num_cores = 8
)

perfil_p <- perfilar_rapido(
    ajuste = aj_iran, df_eventos = df_iran,
    M0 = 4.5,
    x_grid = x_grid_iran, y_grid = y_grid_iran,
    parametro = "p", grade_valores = grade_p,
    num_cores = 8
)

perfil_D <- perfilar_rapido(
    ajuste = aj_iran, df_eventos = df_iran,
    M0 = 4.5,
    x_grid = x_grid_iran, y_grid = y_grid_iran,
    parametro = "D", grade_valores = grade_D,
    num_cores = 8
)

perfil_q <- perfilar_rapido(
    ajuste = aj_iran, df_eventos = df_iran,
    M0 = 4.5,
    x_grid = x_grid_iran, y_grid = y_grid_iran,
    parametro = "q", grade_valores = grade_q,
    num_cores = 8
)

perfil_gamma <- perfilar_rapido(
    ajuste = aj_iran, df_eventos = df_iran,
    M0 = 4.5,
    x_grid = x_grid_iran, y_grid = y_grid_iran,
    parametro = "gamma", grade_valores = grade_gamma,
    num_cores = 8
)

# Perfis:
par(mfrow = c(2, 4))

plot(perfil_A$valor, perfil_A$log_verossimilhanca, type = "l",
    xlab = "A", ylab = "Log-verossimilhança perfilada")
abline(v = aj_iran$theta["A"], lty = 2)

plot(perfil_alpha$valor, perfil_alpha$log_verossimilhanca, type = "l",
    xlab = "alpha", ylab = "Log-verossimilhança perfilada")
abline(v = aj_iran$theta["alpha"], lty = 2)

plot(perfil_c$valor, perfil_c$log_verossimilhanca, type = "l",
    xlab = "c", ylab = "Log-verossimilhança perfilada")
abline(v = aj_iran$theta["c"], lty = 2)

plot(perfil_p$valor, perfil_p$log_verossimilhanca, type = "l",
    xlab = "p", ylab = "Log-verossimilhança perfilada")
abline(v = aj_iran$theta["p"], lty = 2)

plot(perfil_D$valor, perfil_D$log_verossimilhanca, type = "l",
    xlab = "D", ylab = "Log-verossimilhança perfilada")
abline(v = aj_iran$theta["D"], lty = 2)

plot(perfil_q$valor, perfil_q$log_verossimilhanca, type = "l",
    xlab = "q", ylab = "Log-verossimilhança perfilada")
abline(v = aj_iran$theta["q"], lty = 2)

plot(perfil_gamma$valor, perfil_gamma$log_verossimilhanca, type = "l",
    xlab = "gamma", ylab = "Log-verossimilhança perfilada")
abline(v = aj_iran$theta["gamma"], lty = 2)

# Verossimilhanca no MLE (para a estatistica de razão de verossimilhancas):
passo_e_mle <- calcular_probabilidades(df_iran, aj_iran$theta, 4.5,
                                       aj_iran$raster, x_grid_iran, y_grid_iran, 10)

logvero_mle <- -log_lik_etas(aj_iran$theta, df_iran, passo_e_mle$probabilidades, passo_e_mle$prob_total,
                             aj_iran$raster, x_grid_iran, y_grid_iran)

# Estatísticas da RV:
perfil_A$estatistica_rv <- 2 * (logvero_mle - perfil_A$log_verossimilhanca)
perfil_alpha$estatistica_rv <- 2 * (logvero_mle - perfil_alpha$log_verossimilhanca)
perfil_c$estatistica_rv <- 2 * (logvero_mle - perfil_c$log_verossimilhanca)
perfil_p$estatistica_rv <- 2 * (logvero_mle - perfil_p$log_verossimilhanca)
perfil_D$estatistica_rv <- 2 * (logvero_mle - perfil_D$log_verossimilhanca)
perfil_q$estatistica_rv <- 2 * (logvero_mle - perfil_q$log_verossimilhanca)
perfil_gamma$estatistica_rv <- 2 * (logvero_mle - perfil_gamma$log_verossimilhanca)

par(mfrow = c(2, 4))
plot(perfil_A$valor, perfil_A$estatistica_rv, type = "l",
    xlab = "A", ylab = "Estatística RV")
abline(v = aj_iran$theta["A"], lty = 2)
plot(perfil_gamma$valor, perfil_gamma$estatistica_rv, type = "l",
    xlab = "gamma", ylab = "Estatística RV")
abline(v = aj_iran$theta["gamma"], lty = 2)

#====Teste====
# Taxa global de background (Eventos/dia em toda a região):
mu_global_iran <- sum(pmax(0, 1 - aj_iran$prob_total)) / max(df_temp_iran$time)

# Tau_j para os tempos dos eventos:
tau_j_iran <- calcular_tau(df_temp_iran$time, df_temp_iran, aj_iran$theta, mu_global_iran, M0 = 4.5)

# Calcular U_j:
delta_tau_iran <- diff(c(0, tau_j_iran))
U_j_iran <- 1 - exp(-delta_tau_iran)

#----Gráficos----
par(mfrow = c(1, 2))

# Gráfico A: tau_j vs j
plot(1:length(tau_j_iran), tau_j_iran, type = "l", col = "blue", lwd = 2,
     xlab = "Número do Evento (j)", ylab = expression(tau[j]),
     main = "Tempos Transformados")
abline(a = 0, b = 1, col = "red", lty = 2)

# Gráfico B: Q-Q Plot de U_j
plot(sort(U_j_iran), qunif(ppoints(length(U_j_iran))), pch = 20, col = "darkgreen",
     xlab = "Empírico", ylab = "Teórico (Uniforme 0-1)",
     main = "Q-Q Plot de U_j")
abline(a = 0, b = 1, col = "red", lty = 2)

ks.test(U_j_iran, punif)

# Número de caixas (bins) de tempo
n_temp_iran <- 100
janelas_tempo_iran <- seq(0, max(df_temp_iran$time), length.out = n_temp_iran + 1)

# Contagem de eventos em cada janela:
obs_temp_iran <- hist(df_temp_iran$time, breaks = janelas_tempo_iran, plot = FALSE)$counts

# Contagens esperadas pela diferença do tau nos limites da janela:
tau_limites_iran <- calcular_tau(janelas_tempo_iran, df_temp_iran, aj_iran$theta, mu_global_iran, M0 = 4.5)
esp_temp_iran <- diff(tau_limites_iran)

# Resíduos de Pearson Temporais:
res_temp_pearson_iran <- (obs_temp_iran - esp_temp_iran) / sqrt(esp_temp_iran)

# Resíduos Temporais:
par(mfrow = c(1, 1))
plot(janelas_tempo_iran[-1], res_temp_pearson_iran, type = "h", lwd = 2, col = "purple",
     xlab = "Tempo (dias)", ylab = "Resíduos de Pearson",
     main = "Resíduos Temporais Suavizados")
abline(h = 0, col = "black")
abline(h = c(-2, 2), col = "red", lty = 2)

#====Espaciais====

# Definir os pontos observados:
pontos_iran <- ppp(df_temp_iran$x, df_temp_iran$y,
                   window = owin(x_grid_iran[c(1, length(x_grid_iran))],
                                 y_grid_iran[c(1, length(y_grid_iran))]),
              marks = pmax(0, 1 - aj_iran$prob_total))

# Criar a densidade espacial esperada pelo Background:
cell_area_iran <- diff(x_grid_iran)[1] * diff(y_grid_iran)[1]
matriz_esperada_iran <- (aj_iran$raster / (sum(aj_iran$raster) * cell_area_iran)) * (mu_global_iran * max(df_temp_iran$time))
# matriz_esperada_iran <- (rasters_iran[[length(rasters_iran)]] / sum(rasters_iran[[length(rasters_iran)]])) * (mu_global_iran * max(df_temp_iran$time))

# Extrair os 7 parâmetros do modelo convergido
A_iran <- aj_iran$theta["A"]
alpha_iran <- aj_iran$theta["alpha"]
D_iran <- aj_iran$theta["D"]
q_par_iran <- aj_iran$theta["q"]
gamma_par_iran <- aj_iran$theta["gamma"]

M0 <- 4.5 # Ajuste para o limiar de magnitude correto do seu catálogo atual

# Sobrepor a nuvem espacial de Pareto (Lei de Potência) de TODOS os eventos
for(i in 1:nrow(df_temp_iran)) {
    xi_iran <- df_temp_iran$x[i]
    yi_iran <- df_temp_iran$y[i]
    mi_iran <- df_temp_iran$mag[i]
    # Produtividade esperada do evento (mantém-se com alpha)
    k_i_iran <- A_iran * exp(alpha_iran * (mi_iran - M0))
    # Fator de escala espacial S_i (agora independente de alpha, usa gamma)
    S_i_iran <- D_iran * exp(gamma_par_iran * (mi_iran - M0))
    # Matriz de distâncias ao quadrado
    dx2_iran <- (x_grid_iran - xi_iran)^2
    dy2_iran <- (y_grid_iran - yi_iran)^2
    dist2_matrix_iran <- outer(dx2_iran, dy2_iran, "+")
    # Kernel de Pareto Bivariado Isotrópico (Fórmula do pacote ETAS)
    # f(x,y) = ((q-1) / (pi * S_i)) * (1 + dist2 / S_i)^(-q)
    term1_iran <- (q_par_iran - 1) / (pi * S_i_iran)
    term2_iran <- (1 + (dist2_matrix_iran / S_i_iran))^(-q_par_iran)
    kernel_esp_iran <- k_i_iran * term1_iran * term2_iran
    matriz_esperada_iran <- matriz_esperada_iran + kernel_esp_iran
}

densidade_esp_iran <- im(t(matriz_esperada_iran), xcol = x_grid_iran, yrow = y_grid_iran)

# Criar a densidade observada, forçando o alinhamento com xy:
densidade_obs_iran <- density(pontos_iran, sigma = bw.diggle(pontos_iran), xy = densidade_esp_iran,
                              weights = pontos_iran$marks)

# Calcular os resíduos de Pearson:
epsilon <- 1e-10
res_spat_pearson_iran <- eval.im((densidade_obs_iran - densidade_esp_iran) / sqrt(densidade_esp_iran + epsilon))

# Visualização:

min_val <- min(res_spat_pearson_iran)
max_val <- max(res_spat_pearson_iran)

max_abs <- max(abs(c(min_val, max_val)))
sym_breaks <- seq(-max_abs, max_abs, length.out = 257)

# 2. Build a dedicated colourmap object mapping the breaks to the Oslo palette
col_map <- colourmap(
    col = hcl.colors(256, palette = "Viridis", rev = F),
  breaks = sym_breaks
)

plot(res_spat_pearson_iran, main = "Resíduos Espaciais Totais de Pearson (Modelo Pareto)", col = col_map)
plot(res_spat_pearson_iran, main = "Resíduos Espaciais Totais de Pearson (Modelo Pareto)", col = topo.colors(256))

#====Gráficos do pacote ETAS====
cat_dummy <- catalog(iran.quakes)
lat.range <- range(cat_dummy$region.poly$lat)
long.range <- range(cat_dummy$region.poly$long)
xy.bnd <- longlat2xy(c(long.range, rev(long.range)), rep(lat.range,
                                                         each = 2), aj_iran_etas$object$region.poly, aj_iran_etas$object$dist.unit)
rv <- diff(range(xy.bnd$x))/diff(range(xy.bnd$y))
dimyx <- round(128 * c(1, rv))
gx <- seq(min(xy.bnd$x), max(xy.bnd$x), length.out = 128)
gy <- seq(min(xy.bnd$y), max(xy.bnd$y), length.out = dimyx[1])
out <- cxxrates(aj_iran_etas$param, aj_iran_etas$object$revents, aj_iran_etas$bwd, aj_iran_etas$object$rtperiod,
                gx, gy, aj_iran_etas$mver)

grid_iran_2 <- expand.grid(x = seq(min(iran.cat$longlat.coord$long), max(iran.cat$longlat.coord$long), length.out = 128),
                           y = seq(min(iran.cat$longlat.coord$lat), max(iran.cat$longlat.coord$lat), length.out = 128))

par(mfrow = c(2, 1), mar = c(4, 4, 2, 5))

plot(df_iran$x, df_iran$y, xlab = "long",
     ylab = "lat", col = "blue", cex = 2 * (df_iran$mag - 4.4)/max(df_iran$mag - 4.4),
     asp = TRUE, axes = FALSE,
     xlim = range(df_iran$x),
     ylim = range(df_iran$y))
maps::map("world", add = TRUE, col = "grey50")
axis(1)
axis(2)
box()
fields::quilt.plot(grid_iran_2$x, grid_iran_2$y, out$bkgd,
           nx = 64, ny = 64, asp = T,
           xlim = range(df_iran$x),
           ylim = range(df_iran$y),
           add.legend = F)
maps::map("world", add = TRUE, col = "grey50")
axis(1)
axis(2)
box()
fields::image.plot(zlim = range(aj_iran$raster, na.rm = TRUE),
                   legend.only = TRUE)

#====Variabilidade da intensidade bootstrap====
mapa_boot <- mapa_incerteza_fundo_bootstrap(backgrounds_iran)
image(x_grid_iran, y_grid_iran, mapa_boot$dp, main = "Desvio-padrão bootstrap do fundo")
image(x_grid_iran, y_grid_iran, mapa_boot$cv, main = "Coeficiente de variação do fundo")
