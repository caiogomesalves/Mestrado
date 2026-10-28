#====Organização do range====

# Carregar o pacote spatstat
library(spatstat)

# 1. Definir a janela espacial de observação (ex: 10x10)
janela <- owin(xrange = c(-5, 5), yrange = c(-5, 5))

# 2. Definir a função de intensidade de background mu(x, y)
# Exemplo usando uma intensidade não-homogênea com pico no centro (5,5)
mu_func <- function(x, y) {
  0.5 * exp(-(x)^2 - (y)^2)
}

# Definir a dimensão da matriz desejada
n <- 64

imagem_mu <- as.im(mu_func, W = janela, dimyx = c(n, n))
matriz_intensidade <- as.matrix(imagem_mu)

# Visualização rápida para o seu pôster via spatstat base:
plot(imagem_mu, main = "Matriz de Intensidade mu(s)")

#====Simulação====

theta_verdadeiro <- c("A" = 0.5, "alpha" = 1.5, "c" = 0.5,
                      "p" = 1.5, "D" = 0.1, "q" = 2, "gamma" = 1)

x_grid_sim <- seq(-5, 5, length.out = 64)
y_grid_sim <- seq(-5, 5, length.out = 64)

dados_simulados_poster <- simular_catalogo_etas(theta = theta_verdadeiro, M0 = 4.5, b = 2, Mmax = 10,
                                                T_max = 365, mu_type = "raster", mu_raster = matriz_intensidade,
                                                x_grid = x_grid_sim, y_grid = y_grid_sim, semente = 0)

cat_simulado_poster <- catalog(
    dados_simulados_poster %>%
    mutate(t = ddays(time), data = POSIXct(1) + t,
           long = x, lat = y, date = date(data),
           time = format(data, format = "%H:%M:%S"))
)

#====Ajustes====

t1_simulado <- Sys.time()

aj_poster <- ajustar_etas_proprio(catalogo = cat_simulado_poster, theta0 = theta_verdadeiro,
                                  M0 = 4.5, x_grid = x_grid_sim, y_grid = y_grid_sim)

t2_simulado <- Sys.time()

t3_simulado <- Sys.time()

aj_poster_etas <- etas(cat_simulado_poster, param0 = c("mu" = 1, theta_verdadeiro))

t4_simulado <- Sys.time()

#====Resíduos====
# Taxa global de background (Eventos/dia em toda a região):
mu_global_simulado <- sum(pmax(0, 1 - aj_poster$prob_total)) / max(cat_simulado_poster$rtperiod)

# Tau_j para os tempos dos eventos:
tau_j_simulado <- calcular_tau(dados_simulados_poster$time, dados_simulados, aj_poster$theta, mu_global_simulado, M0 = 4.5)

# Calcular U_j:
delta_tau_simulado <- diff(c(0, tau_j_simulado))
U_j_simulado <- 1 - exp(-delta_tau_simulado)

#----Gráficos----
par(mfrow = c(1, 2))

# Gráfico A: tau_j vs j
plot(1:length(tau_j_simulado), tau_j_simulado, type = "l", col = "blue", lwd = 2,
     xlab = "Número do Evento (j)", ylab = expression(tau[j]),
     main = "Tempos Transformados")
abline(a = 0, b = 1, col = "red", lty = 2)

# Gráfico B: Q-Q Plot de U_j
plot(sort(U_j_simulado), qunif(ppoints(length(U_j_simulado))), pch = 20, col = "darkgreen",
     xlab = "Empírico", ylab = "Teórico (Uniforme 0-1)",
     main = "Q-Q Plot de U_j")
abline(a = 0, b = 1, col = "red", lty = 2)

# Número de caixas (bins) de tempo
n_temp <- 100
janelas_tempo_simulado <- seq(0, max(dados_simulados_poster$time), length.out = n_temp + 1)

# Contagem de eventos em cada janela:
obs_temp_simulado <- hist(dados_simulados_poster$time, breaks = janelas_tempo_simulado, plot = FALSE)$counts

# Contagens esperadas pela diferença do tau nos limites da janela:
tau_limites_simulado <- calcular_tau(janelas_tempo_simulado, dados_simulados_poster, aj_poster$theta, mu_global_simulado, M0 = 4.5)
esp_temp_simulado <- diff(tau_limites_simulado)

# Resíduos de Pearson Temporais:
res_temp_pearson_simulado <- (obs_temp_simulado - esp_temp_simulado) / sqrt(esp_temp_simulado)

# Resíduos Temporais:
par(mfrow = c(1, 1))
plot(janelas_tempo_simulado[-1], res_temp_pearson_simulado, type = "h", lwd = 2, col = "purple",
     xlab = "Tempo (dias)", ylab = "Resíduos de Pearson",
     main = "Resíduos Temporais Suavizados")
abline(h = 0, col = "black")
abline(h = c(-2, 2), col = "red", lty = 2)

# Definir os pontos observados:
pontos_simulado <- ppp(dados_simulados_poster$x, dados_simulados_poster$y,
                       window = owin(x_grid_sim[c(1, length(x_grid_sim))], y_grid_sim[c(1, length(y_grid_sim))]),
                       marks = pmax(0, 1 - aj_poster$prob_total))

# Criar a densidade espacial esperada:
# matriz_esperada <- (rasters[[length(rasters)]] / sum(rasters[[length(rasters)]])) * (mu_global * max(df_temp$time)
cell_area_simulado <- diff(x_grid_sim)[1] * diff(y_grid_sim)[1]
matriz_esperada_sim <- (aj_poster$raster / (sum(aj_poster$raster) * cell_area_simulado)) * (mu_global_simulado * max(dados_simulados_poster$time))
# matriz_esperada <- rasters[[length(rasters)]]

A_sim <- aj_poster$theta["A"]
alpha_sim <- aj_poster$theta["alpha"]
D_sim <- aj_poster$theta["D"]
q_par_sim <- aj_poster$theta["q"]
gamma_par_sim <- aj_poster$theta["gamma"]

M0 <- 4.5 # Ajuste para o limiar de magnitude correto do seu catálogo atual

# Sobrepor a nuvem espacial de Pareto (Lei de Potência) de TODOS os eventos
for(i in 1:nrow(dados_simulados_poster)) {
    xi_sim <- dados_simulados_poster$x[i]
    yi_sim <- dados_simulados_poster$y[i]
    mi_sim <- dados_simulados_poster$mag[i]
    # Produtividade esperada do evento (mantém-se com alpha)
    k_i_sim <- A_sim * exp(alpha_sim * (mi_sim - M0))
    # Fator de escala espacial S_i (agora independente de alpha, usa gamma)
    S_i_sim <- D_sim * exp(gamma_par_sim * (mi_sim - M0))
    # Matriz de distâncias ao quadrado
    dx2_sim <- (x_grid_sim - xi_sim)^2
    dy2_sim <- (y_grid_sim - yi_sim)^2
    dist2_matrix_sim <- outer(dx2_sim, dy2_sim, "+")
    # Kernel de Pareto Bivariado Isotrópico (Fórmula do pacote ETAS)
    # f(x,y) = ((q-1) / (pi * S_i)) * (1 + dist2 / S_i)^(-q)
    term1_sim <- (q_par_sim - 1) / (pi * S_i_sim)
    term2_sim <- (1 + (dist2_matrix_sim / S_i_sim))^(-q_par_sim)
    kernel_esp_sim <- k_i_sim * term1_sim * term2_sim
    matriz_esperada_sim <- matriz_esperada_sim + kernel_esp_sim
}

densidade_esp_sim <- im(t(matriz_esperada_sim), xcol = x_grid_sim, yrow = y_grid_sim)

# Criar a densidade observada, forçando o alinhamento com xy:
densidade_obs_sim <- density(pontos_simulado, sigma = bw.diggle(pontos_simulado), xy = densidade_esp_sim,
                             weights = pontos_simulado$marks)

# Calcular os resíduos de Pearson:
epsilon <- 1e-10
res_spat_pearson_sim <- eval.im((densidade_obs_sim - densidade_esp_sim) / sqrt(densidade_esp_sim + epsilon))

# Visualização:
plot(res_spat_pearson_sim, main = "Resíduos Espaciais Totais de Pearson (Modelo Pareto)", col = topo.colors(256))
contour(res_spat_pearson_sim, add = TRUE, col = "black", alpha = 0.5)
