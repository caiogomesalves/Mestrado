#----Funções auxiliares----

# Produtividade:
kappa_func <- function(m, m0, A, alpha) {
    return(A * exp(alpha * (m - m0)))
}

# Densidade de probabilidade do tempo de ocorrência (Lei de Omori modificada):
g_func <- function(t_diff, c, p) {
    val <- rep(0, length(t_diff))
    idx <- t_diff > 0
    val[idx] <- (p - 1) * (c^(p - 1)) * (t_diff[idx] + c)^(-p)
    return(val)
}

# Densidade de probabilidade da localização espacial:
f_func <- function(x_diff, y_diff, m, m0, D, gamma, q) {
    # Distância ao quadrado
    r2 <- x_diff^2 + y_diff^2
    escala <- D * exp(gamma * (m - m0))
    densidade <- ((q - 1) / (pi * escala)) * (1 + r2 / escala)^(-q)
    return(densidade)
}

# Integral da distribuição temporal g(t) no intervalo [t_start, t_start + T]:
integral_g_func <- function(t_i, t_start, T_len, c, p) {
    # Conforme equação 6 e 130 do artigo
    int_val <- numeric(length(t_i))
    idx_before <- t_i <= t_start
    idx_after <- t_i > t_start
    # Para eventos que ocorreram antes do início do período de estudo:
    if(any(idx_before)) {
        t_i_b <- t_i[idx_before]
        int_val[idx_before] <- (1 + (t_start - t_i_b)/c)^(1 - p) -
            (1 + (t_start + T_len - t_i_b)/c)^(1 - p)
    }
    # Para eventos que ocorreram durante o período de estudo:
    if(any(idx_after)) {
        t_i_a <- t_i[idx_after]
        int_val[idx_after] <- 1 - (1 + (t_start + T_len - t_i_a)/c)^(1 - p)
    }
    return(int_val)
}

#----LogVerossimilhança L1----

l1_loglik <- function(beta, cat_obj) {
    # Extraindo dados do objeto catalog:
    m0 <- cat_obj$mag.threshold
    m <- cat_obj$revents[, "mm"] + m0

    # Determinando eventos alvo, com flag == 1:
    delta <- cat_obj$revents[, "flag"]

    N_prime <- sum(delta)

    # Equação 5 (l1)
    l1 <- log(beta) * N_prime - beta * sum(delta * (m - m0))
    return(l1)
}

#----LogVerossimilhança L2----

l2_loglik <- function(theta, cat_obj, u_vals = NULL) {
    # theta = (mu, A, alpha, c, p, D, gamma, q)
    mu <- theta[1]; A <- theta[2]; alpha <- theta[3]; c <- theta[4]
    p <- theta[5]; D <- theta[6]; gamma <- theta[7]; q <- theta[8]

    # Extraindo informações do catálogo:
    t <- cat_obj$revents[, "tt"]
    x <- cat_obj$revents[, "xx"]
    y <- cat_obj$revents[, "yy"]
    m0 <- cat_obj$mag.threshold
    m <- cat_obj$revents[, "mm"] + m0

    # Tempos de estudo (em dias decorridos):
    t_start <- t[1]
    T_len <- t[length(t)]

    N <- length(t)

    # Determinando eventos alvo (delta_i = 1) e complementares (delta_i = 0):
    delta <- cat_obj$revents[, "flag"]

    # Se o background espacial u(x,y) não for fornecido, assumimos uniforme = 1/Área
    if (is.null(u_vals)) {
        u_vals <- rep(1, N)
        integral_u <- 1
    } else {
        integral_u <- 1
    }

    # 1. CÁLCULO DO TERMO DE SOMATÓRIO (para eventos alvo)
    sum_log_lambda <- 0

    for (j in 1:N) {
        if (delta[j] == 1) {
            # Eventos que ocorreram antes do evento j
            hist_idx <- t < t[j]

            if (any(hist_idx)) {
                t_hist <- t[hist_idx]
                x_hist <- x[hist_idx]
                y_hist <- y[hist_idx]
                m_hist <- m[hist_idx]

                # Componentes da intensidade:
                kappa_val <- kappa_func(m_hist, m0, A, alpha)
                g_val <- g_func(t[j] - t_hist, c, p)
                f_val <- f_func(x[j] - x_hist, y[j] - y_hist, m_hist, m0, D, gamma, q)

                # Intensidade condicional no ponto j:
                lambda_j <- mu * u_vals[j] + sum(kappa_val * g_val * f_val)

                sum_log_lambda <- sum_log_lambda + log(lambda_j)
            } else {
                # Se não há histórico, a intensidade é apenas o background
                sum_log_lambda <- sum_log_lambda + log(mu * u_vals[j])
            }
        }
    }

    # 2. CÁLCULO DO TERMO DE INTEGRAL (sobre o espaço e o tempo)
    # Integral do background
    int_background <- T_len * mu * integral_u

    # Integral das réplicas
    int_g <- integral_g_func(t, t_start, T_len, c, p)
    kappa_all <- kappa_func(m, m0, A, alpha)

    # Aproximação: assumindo que a integral espacial de f(x,y) sobre a região S é ~1:
    int_f <- rep(1, N)

    int_trigger <- sum(kappa_all * int_g * int_f)

    # Equação 5 (l2)
    l2 <- sum_log_lambda - (int_background + int_trigger)

    # Retornamos o log-likelihood negativo para usar funções de minimização como optim()
    return(-l2)
}

#----Aplicação nos dados simulados----

##----Beta----

# Dados simulados:
plot(simu_peru)

# Cálculo analítico do MLE para beta:
delta <- rep(1, nrow(simu_peru$revents))

N_prime <- sum(delta)

beta_est <- N_prime / sum(delta * (simu_peru$revents[, "mm"]))

# Estimativa por maximização da função de verossimilhança:
beta_est_loglik <- optim(par = c(beta = 1), function(x){-l1_loglik(x, simu_peru)},
                         method = "Brent",
                         lower = 0.1, upper = 10)$par

curve(-l1_loglik(x, simu_peru), 0.1, 10,
      xlab = expression(mu))

##----Theta----

# Mu:
optim(par = c(mu = 0.15, theta_peru), fn = l2_loglik,
      cat_obj = simu_peru,
      method = "L-BFGS-B",
      lower = c(0.1, 0.05, 0.05, 1.01, 0.001, 0.001, 1.01, 2), upper = 5)

mu_est <- apply(matrix(seq(0.1, 2, length.out = 100), ncol = 1), 1,
                function(x){l2_loglik(c(x, theta_peru), simu_peru)})

# Gráfico da verossimilhança:
plot(x = seq(0.1, 2, length.out = 100),
     mu_est,
     type = "l", xlab = expression(mu), ylab = "Log-Verossimilhança")

# Outros parâmetros:
optim(par = c(A = 1), function(x) {l2_loglik(c(mu = 0.9582442, x, theta_peru[-1]), simu_peru)},
      method = "Brent",
      lower = 0.01, upper = 5)

seq_x <- seq(-11, 11, length.out = 100)
val_l2 <- numeric(100)

for (i in 1:length(seq_x)) {
    val_l2[i] <- l2_loglik(c(0.9582442, seq_x[i], theta_peru[-1]), simu_peru)
}

plot(seq_x, val_l2, type = "l")

# Dois a dois:
grid <- expand.grid(mu = seq(0.01, 5, length.out = 10),
                    A = seq(0.01, 5, length.out = 10))

loglik_teste <- apply(grid, 1, function(x) {l2_loglik(c(x[1], x[2],
                                                        theta_peru[-c(1)]),
                                                      simu_peru)})

grid$loglik <- loglik_teste

ggplot(grid, aes(x = mu, y = A, fill = loglik)) +
    geom_tile() +
    coord_equal()
