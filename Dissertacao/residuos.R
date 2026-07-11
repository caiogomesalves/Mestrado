#----Pacotes----

library(tidyverse)
library(ETAS)
library(spatstat)
library(readxl)
library(parallel)
library(terra)

#----Resultados----

load("~/resultados_thinning_final.RData")

#----Funções----

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

#====Tempos Transformados====

# Taxa global de background (Eventos/dia em toda a região):
mu_global <- sum(pmax(0, 1 - prob_total)) / max(df_temp$time)

# Tau_j para os tempos dos eventos:
tau_j <- calcular_tau(df_temp$time, df_temp, thetas[[length(thetas)]], mu_global, M0 = 4.5)

# Calcular U_j:
delta_tau <- diff(c(0, tau_j))
U_j <- 1 - exp(-delta_tau)

#----Gráficos----
par(mfrow = c(1, 2))

# Gráfico A: tau_j vs j
plot(1:length(tau_j), tau_j, type = "l", col = "blue", lwd = 2,
     xlab = "Número do Evento (j)", ylab = expression(tau[j]),
     main = "Tempos Transformados")
abline(a = 0, b = 1, col = "red", lty = 2)

# Gráfico B: Q-Q Plot de U_j
plot(sort(U_j), qunif(ppoints(length(U_j))), pch = 20, col = "darkgreen",
     xlab = "Empírico", ylab = "Teórico (Uniforme 0-1)",
     main = "Q-Q Plot de U_j")
abline(a = 0, b = 1, col = "red", lty = 2)

#====Temporais====

# Número de caixas (bins) de tempo
n_temp <- 100
janelas_tempo <- seq(0, max(df_temp$time), length.out = n_temp + 1)

# Contagem de eventos em cada janela:
obs_temp <- hist(df_temp$time, breaks = janelas_tempo, plot = FALSE)$counts

# Contagens esperadas pela diferença do tau nos limites da janela:
tau_limites <- calcular_tau(janelas_tempo, df_temp, thetas[[length(thetas)]], mu_global, M0 = 4.5)
esp_temp <- diff(tau_limites)

# Resíduos de Pearson Temporais:
res_temp_pearson <- (obs_temp - esp_temp) / sqrt(esp_temp)

# Resíduos Temporais:
par(mfrow = c(1, 1))
plot(janelas_tempo[-1], res_temp_pearson, type = "h", lwd = 2, col = "purple",
     xlab = "Tempo (dias)", ylab = "Resíduos de Pearson",
     main = "Resíduos Temporais Suavizados")
abline(h = 0, col = "black")
abline(h = c(-2, 2), col = "red", lty = 2)

#====Espaciais====

# Definir os pontos observados:
pontos <- ppp(df_temp$x, df_temp$y,
              window = owin(x_grid[c(1, length(x_grid))], y_grid[c(1, length(y_grid))]),
              marks = pmax(0, 1 - prob_total))

# Criar a densidade espacial esperada:
# matriz_esperada <- (rasters[[length(rasters)]] / sum(rasters[[length(rasters)]])) * (mu_global * max(df_temp$time)
cell_area <- diff(x_grid)[1] * diff(y_grid)[1]
matriz_esperada <- (rasters[[length(rasters)]] / (sum(rasters[[length(rasters)]]) * cell_area)) * (mu_global * max(df_temp$time))
# matriz_esperada <- rasters[[length(rasters)]]

A <- thetas[[length(thetas)]]["A"]
alpha <- thetas[[length(thetas)]]["alpha"]
D <- thetas[[length(thetas)]]["D"]
q_par <- thetas[[length(thetas)]]["q"]
gamma_par <- thetas[[length(thetas)]]["gamma"]

M0 <- 4.5 # Ajuste para o limiar de magnitude correto do seu catálogo atual

# Sobrepor a nuvem espacial de Pareto (Lei de Potência) de TODOS os eventos
for(i in 1:nrow(df_temp)) {
    xi <- df_temp$x[i]
    yi <- df_temp$y[i]
    mi <- df_temp$mag[i]
    # Produtividade esperada do evento (mantém-se com alpha)
    k_i <- A * exp(alpha * (mi - M0))
    # Fator de escala espacial S_i (agora independente de alpha, usa gamma)
    S_i <- D * exp(gamma_par * (mi - M0))
    # Matriz de distâncias ao quadrado
    dx2 <- (x_grid - xi)^2
    dy2 <- (y_grid - yi)^2
    dist2_matrix <- outer(dx2, dy2, "+")
    # Kernel de Pareto Bivariado Isotrópico (Fórmula do pacote ETAS)
    # f(x,y) = ((q-1) / (pi * S_i)) * (1 + dist2 / S_i)^(-q)
    term1 <- (q_par - 1) / (pi * S_i)
    term2 <- (1 + (dist2_matrix / S_i))^(-q_par)
    kernel_esp <- k_i * term1 * term2
    matriz_esperada <- matriz_esperada + kernel_esp
}

densidade_esp <- im(t(matriz_esperada), xcol = x_grid, yrow = y_grid)

# Criar a densidade observada, forçando o alinhamento com xy:
densidade_obs <- density(pontos, sigma = bw.diggle(pontos), xy = densidade_esp,
                         weights = pontos$marks)

# Calcular os resíduos de Pearson:
epsilon <- 1e-10
res_spat_pearson <- eval.im((densidade_obs - densidade_esp) / sqrt(densidade_esp + epsilon))

# Visualização:
plot(res_spat_pearson, main = "Resíduos Espaciais Totais de Pearson (Modelo Pareto)", col = topo.colors(256))
contour(res_spat_pearson, add = TRUE, col = "black", alpha = 0.5)

#====Tolerância====

diferencas_rasters <- list()

for (i in 2:length(rasters)) {
    diferencas_rasters[[i]] <- max(abs(rasters[[i]] - rasters[[i - 1]]))
}

diferencas_thetas <- list()

for (i in 2:length(thetas)) {
    diferencas_thetas[[i]] <- max(abs(thetas[[i]] - thetas[[i - 1]]))
}

do.call(c, diferencas_rasters) %>%
    cbind(do.call(c, diferencas_thetas)) %>%
    matplot(type = "b", pch = c(1, 1))
legend("topright", legend = c("Rasteres", "Parâmetros"),
       col = c("black", "red"), pch = c(1, 1))
abline(h = 0.01)

ts(do.call(rbind, thetas)) %>%
    lattice::xyplot(xlab = "Iteração", pch = 16, type = "b",
                    scales = list(x = list(tick.number = length(thetas)),
                                  y = list(tick.number = 3)),
                    main = "Estimativa dos parâmetros do modelo")

#====Resíduos ponderados====

# Definir os pontos com os pesos (Probabilidade de ser Background):
pesos_bg <- pmax(0, 1 - prob_total)

pontos_bg <- ppp(df_temp$x, df_temp$y,
                 window = owin(x_grid[c(1, length(x_grid))], y_grid[c(1, length(y_grid))]),
                 marks = pesos_bg)

# Densidade Esperada é o raster base normalizado:
cell_area <- diff(x_grid)[1] * diff(y_grid)[1]
matriz_esperada_bg <- (rasters[[length(rasters)]] / sum(rasters[[length(rasters)]]) * cell_area) * (mu_global * max(df_temp$time))
densidade_esp_bg <- im(t(matriz_esperada_bg), xcol = x_grid, yrow = y_grid)

# Densidade observada do background:
densidade_obs_bg <- density(pontos_bg,
                            weights = pontos_bg$marks,
                            sigma = bw.diggle(pontos_bg),
                            xy = densidade_esp_bg)

# Calcular os resíduos de Pearson do background:
epsilon <- 1e-10
res_bg_pearson <- eval.im((densidade_obs_bg - densidade_esp_bg) / sqrt(densidade_esp_bg + epsilon))

# Visualização:
plot(res_bg_pearson, main = "Resíduos de Pearson (Apenas Background)", col = topo.colors(256))
contour(res_bg_pearson, add = TRUE, col = "black", alpha = 0.5)


#====Erros-padrão====

hessianas[[length(hessianas)]]

matriz_covariancia <- solve(hessianas[[length(hessianas)]])

erros_padrao <- sqrt(diag(matriz_covariancia))

#====Simulação com pacote====

bkgd <- list(x = x_grid,
             y = y_grid,
             bkgd = rasters[[length(rasters)]])

simetas(c(thetas[[length(thetas)]], "beta" = 3.122028),
        bkgd,
        sim.start = 0,
        sim.length = 10000,
        mag.threshold = 4.5,
        region.poly = peru.quakes$region.poly) %>%
    plot()
