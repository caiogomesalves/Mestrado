#----Pacotes----

library(ETAS)

library(tidyverse)

#----Simulador----

simulate_hawkes_etas_jss <- function(
                                     # Parâmetros gerais:
                                     T_max = 100,           # Tempo máximo
                                     X_lim = c(-10, 10),    # Limites em X
                                     Y_lim = c(-10, 10),    # Limites em Y
                                     mu = 0.5,              # Taxa de background base
                                     u_fn = NULL,           # Função não-homogênea u(x,y) que retorna a intensidade espacial
                                     u_max = 1.0,           # Valor máximo da função u_fn na janela para aceitação/rejeição
                                     # Parâmetros ETAS:
                                     A = 0.1,               # Produtividade
                                     alpha = 1.05,          # Sensibilidade da produtividade à magnitude
                                     c_par = 0.05,          # Deslocamento temporal
                                     p_par = 1.2,           # Decaimento temporal (p > 1)
                                     D_par = 0.02,          # Escala espacial base
                                     gamma = 0.5,           # Efeito da magnitude na escala espacial
                                     q_par = 1.5,           # Decaimento espacial (q > 1)
                                     beta = 1.5,            # Parâmetro da exponencial para distribuição das magnitudes
                                     m_min = 2.0            # Magnitude mínima (m0)
                                     ){
    # Lista para armazenar todos os eventos:
    all_events <- data.frame(t = numeric(), x = numeric(), y = numeric(), m = numeric(), gen = numeric())
    area <- (X_lim[2] - X_lim[1]) * (Y_lim[2] - Y_lim[1])
    # ---- 1. Geração dos Eventos Pais (Background Não-Homogêneo) ----
    # Algoritmo de Thinning (Aceitação-Rejeição)
    # Caso seja especificada função de intensidade:
    if (!is.null(u_fn)) {
        lambda_max <- mu * u_max * area * T_max
        n_candidates <- rpois(1, lambda = lambda_max)
        if (n_candidates > 0) {
            c_x <- runif(n_candidates, X_lim[1], X_lim[2])
            c_y <- runif(n_candidates, Y_lim[1], Y_lim[2])
            mat_pontos <- matrix(c(c_x, c_y), ncol = 2)
            # Probabilidade de aceitar o evento baseada em u(x,y)
            prob_accept <- terra::extract(u_fn, mat_pontos) / u_max
            accepted <- runif(n_candidates) < prob_accept
            n_parents <- sum(accepted)
            if (n_parents > 0) {
                parents <- data.frame(
                    t = runif(n_parents, 0, T_max),
                    x = c_x[accepted],
                    y = c_y[accepted],
                    m = rexp(n_parents, rate = beta) + m_min,
                    gen = 0,
                    parent = 1:n_parents
                )
                all_events <- rbind(all_events, parents)
                current_gen <- parents
            } else {
                return(NULL)
            }
        } else {
            return(NULL)
        }
    } else {
        # Processo de Poisson Homogêneo se u_fn não for fornecido
        n_parents <- rpois(1, lambda = mu * area * T_max)
        if (n_parents > 0) {
            parents <- data.frame(
                t = runif(n_parents, 0, T_max),
                x = runif(n_parents, X_lim[1], X_lim[2]),
                y = runif(n_parents, Y_lim[1], Y_lim[2]),
                m = rexp(n_parents, rate = beta) + m_min,
                gen = 0,
                parent = 1:n_parents
            )
            all_events <- rbind(all_events, parents)
            current_gen <- parents
        } else {
            return(NULL)
        }
    }
    # ---- 2. Geração de Descendentes ----
    generation <- 0
    while(nrow(current_gen) > 0) {
        generation <- generation + 1
        next_gen <- data.frame()
        for(i in 1:nrow(current_gen)) {
            parent <- current_gen[i,]
            # Produtividade: N ~ Poisson( A * exp(alpha * (m - m_min)) )
            expected_offspring <- A * exp(alpha * (parent$m - m_min))
            n_offspring <- rpois(1, expected_offspring)
            if(n_offspring > 0) {
                # 2a. Geração de Tempos (Inversa da CDF de Omori modificada)
                U_t <- runif(n_offspring)
                dt <- c_par * ((1 - U_t)^(-1 / (p_par - 1)) - 1)
                t_new <- parent$t + dt
                valid_t <- t_new <= T_max
                n_val <- sum(valid_t)
                if(n_val > 0) {
                    t_new <- t_new[valid_t]
                    # 2b. Geração do Espaço (Inversa da CDF da distribuição pareto)
                    # Variância espacial dependente da magnitude
                    S_val <- D_par * exp(gamma * (parent$m - m_min))
                    U_s <- runif(n_val)
                    R_squared <- S_val * ((1 - U_s)^(-1 / (q_par - 1)) - 1)
                    R <- sqrt(R_squared)
                    theta <- runif(n_val, 0, 2 * pi)
                    x_new <- parent$x + R * cos(theta)
                    y_new <- parent$y + R * sin(theta)
                    # 2c. Marcas / Magnitudes
                    m_new <- rexp(n_val, rate = beta) + m_min
                    offspring <- data.frame(
                        t = t_new,
                        x = x_new,
                        y = y_new,
                        m = m_new,
                        gen = generation,
                        parent = i
                    )
                    if (!is.null(u_fn)) {
                        prob_accept_n <- terra::extract(u_fn, offspring[, c("x", "y")])$lyr.1 / u_max
                        prob_accept_n[is.na(prob_accept_n)] <- 0
                        accepted_n <- runif(nrow(offspring)) < prob_accept_n
                        offspring <- offspring[accepted_n, ]
                    }
                    # Filtro espacial para manter na janela de observação
                    in_window <- offspring$x >= X_lim[1] & offspring$x <= X_lim[2] &
                        offspring$y >= Y_lim[1] & offspring$y <= Y_lim[2]
                    if(sum(in_window, na.rm = T) > 0) {
                        next_gen <- rbind(next_gen, offspring[in_window,])
                    }
                }
            }
        }
        if(nrow(next_gen) > 0) {
            all_events <- rbind(all_events, next_gen)
            current_gen <- next_gen
        } else {
            break
        }
        # Prevenção de loops infinitos para regimes explosivos/supercríticos
        if(generation > 100){
            warning("Regime supercrítico detectado (mais de 100 gerações). Interrompendo.")
            break
        }
    }
    return(all_events[order(all_events$t), ])
}

simulate_hawkes_etas_zhuang <- function(
                                        # Parâmetros gerais:
                                        T_max = 100,           # Tempo máximo
                                        X_lim = c(-10, 10),    # Limites em X
                                        Y_lim = c(-10, 10),    # Limites em Y
                                        mu = 0.5,              # Taxa de background base
                                        u_fn = NULL,           # Função u(x,y) para o background não-homogêneo
                                        u_max = 1.0,           # Valor máximo de u_fn para aceitação/rejeição
                                        # Parâmetros ETAS (Zhuang 2002):
                                        A = 0.1,               # Parâmetro de produtividade de escala
                                        alpha = 1.05,          # Eficiência de magnitude (Tempo e Espaço)
                                        c_par = 0.05,          # Deslocamento temporal (c)
                                        p_par = 1.2,           # Decaimento temporal de Omori (p > 1)
                                        d_par = 0.02,          # Escala de dispersão espacial (d)
                                        beta = 1.5,            # Parâmetro do Gutenberg-Richter
                                        m_min = 2.0            # Magnitude mínima (M0)
                                        ){
    # Lista para armazenar todos os eventos:
    all_events <- data.frame(t = numeric(), x = numeric(), y = numeric(), m = numeric(), gen = numeric())
    area <- (X_lim[2] - X_lim[1]) * (Y_lim[2] - Y_lim[1])
    # ---- 1. Geração dos Eventos Pais (Background / Imigrantes) ----
    if (!is.null(u_fn)) {
        lambda_max <- mu * u_max * area * T_max
        n_candidates <- rpois(1, lambda = lambda_max)
        if (n_candidates > 0) {
            c_x <- runif(n_candidates, X_lim[1], X_lim[2])
            c_y <- runif(n_candidates, Y_lim[1], Y_lim[2])
            mat_pontos <- matrix(c(c_x, c_y), ncol = 2)
            prob_accept <- terra::extract(u_fn, mat_pontos) / u_max
            prob_accept[is.na(prob_accept)] <- 0
            accepted <- runif(n_candidates) < prob_accept
            n_parents <- sum(accepted)
            if (n_parents > 0) {
                parents <- data.frame(
                    t = runif(n_parents, 0, T_max),
                    x = c_x[accepted],
                    y = c_y[accepted],
                    m = rexp(n_parents, rate = beta) + m_min,
                    gen = 0,
                    parent = 1:n_parents
                )
                all_events <- rbind(all_events, parents)
                current_gen <- parents
            } else {
                return(NULL)
            }
        } else {
            return(NULL)
        }
    } else {
        n_parents <- rpois(1, lambda = mu * area * T_max)
        if (n_parents > 0) {
            parents <- data.frame(
                t = runif(n_parents, 0, T_max),
                x = runif(n_parents, X_lim[1], X_lim[2]),
                y = runif(n_parents, Y_lim[1], Y_lim[2]),
                m = rexp(n_parents, rate = beta) + m_min,
                gen = 0,
                parent = 1:n_parents
            )
            all_events <- rbind(all_events, parents)
            current_gen <- parents
        } else {
            return(NULL)
        }
    }
    # ---- 2. Geração de Descendentes (Mecanismo Branching) ----
    generation <- 0
    while(nrow(current_gen) > 0) {
        generation <- generation + 1
        next_gen <- data.frame()
        for(i in 1:nrow(current_gen)) {
            parent <- current_gen[i,]
            # Produtividade: kappa(M) = A * exp(alpha * (M - M0))
            expected_offspring <- A * exp(alpha * (parent$m - m_min))
            n_offspring <- rpois(1, expected_offspring)
            if(n_offspring > 0) {
                # 2a. Geração de Tempos (Inversa da CDF de Omori-Utsu)
                U_t <- runif(n_offspring)
                dt <- c_par * ((1 - U_t)^(-1 / (p_par - 1)) - 1)
                t_new <- parent$t + dt
                # Filtrar imediatamente os que passam do tempo máximo histórico
                valid_t <- t_new <= T_max
                n_val <- sum(valid_t)
                if(n_val > 0) {
                    t_new <- t_new[valid_t]
                    # 2b. Geração do Espaço (Distribuição Gaussiana de Zhuang 2002)
                    # Variância sigma^2 depende de d_par e exponencialmente de alpha e da magnitude
                    sigma2 <- d_par * exp(alpha * (parent$m - m_min))
                    sd_spatial <- sqrt(sigma2)
                    # Sorteio independente em X e Y simulando uma normal bivariada isotrópica
                    x_new <- rnorm(n_val, mean = parent$x, sd = sd_spatial)
                    y_new <- rnorm(n_val, mean = parent$y, sd = sd_spatial)
                    # 2c. Marcas / Magnitudes (Assunção e: Independentes do pai)
                    m_new <- rexp(n_val, rate = beta) + m_min
                    offspring <- data.frame(
                        t = t_new,
                        x = x_new,
                        y = y_new,
                        m = m_new,
                        gen = generation,
                        parent = parent$parent
                    )
                    # Filtro espacial para manter estritamente dentro da janela de interesse
                    in_window <- offspring$x >= X_lim[1] & offspring$x <= X_lim[2] &
                                 offspring$y >= Y_lim[1] & offspring$y <= Y_lim[2]
                    if(sum(in_window) > 0) {
                        next_gen <- rbind(next_gen, offspring[in_window,])
                    }
                }
            }
        }
        if(nrow(next_gen) > 0) {
            all_events <- rbind(all_events, next_gen)
            current_gen <- next_gen
        } else {
            break
        }
        if(generation > 100){
            warning("Regime supercrítico detectado (mais de 100 gerações). Interrompendo simulação.")
            break
        }
    }
    if(nrow(all_events) == 0) return(NULL)
    return(all_events[order(all_events$t), ])
}

# Definindo uma função de background homogênea:
u_fn_homog <- function(x, y) {
    return(1)
}

# Definindo uma função de background não-homogênea:
u_fn_inhomog <- function(x, y) {
  return(exp(-0.05 * (x^2 + y^2)))
}

val_x <- seq(-10, 10, length.out = 100)
val_y <- seq(-10, 10, length.out = 100)

grid <- expand.grid(val_x, val_y)

val_z_homog <- apply(grid, 1, function(x){u_fn_homog(x[1], x[2])})
val_z_inhomog <- apply(grid, 1, function(x){u_fn_inhomog(x[1], x[2])})

fields::image.plot(val_x,
                   val_y,
                   matrix(val_z_homog,
                          nrow = 100))

fields::image.plot(val_x,
                   val_y,
                   matrix(val_z_inhomog,
                          nrow = 100))

#----Simulação----

set.seed(5050)

# Homogêneo:

sim_dados_homog <- simulate_hawkes_etas_jss(
  T_max = 20,
  X_lim = c(-1, 1), Y_lim = c(-1, 1),
  mu = 10, u_fn = u_fn_homog, u_max = 1,
  A = 2.1, alpha = 2.1, c_par = 1.5, p_par = 2,
  D_par = 0.5, gamma = 0.5, q_par = 1, beta = 2, m_min = 4
)

catalogo_dados_homog <- catalog(
    sim_dados_homog %>%
    mutate(t = dyears(t),
           data = POSIXct(1) + t,
           long = x, lat = y,
           mag = m,
           date = date(data),
           time = format(data, format = "%H:%M:%S"))
)

plot(catalogo_dados_homog)

# Não roda:
aj_teste_homog_1 <- etas(catalogo_dados_homog, nthreads = 8)

aj_teste_homog_2 <- etas(catalogo_teste,
                         param0 = c(mu = 0.05, A = 1.1,
                                    alpha = 1.1, c = 1.5,
                                    p = 2, D = 0.5, gamma = 0.5,
                                    q = 1.15),
                         nthreads = 8)

ETAS::simetas()

plot(aj_teste_homog_1)

rates(aj_teste_homog_1)

# Não-homogêneo:

set.seed(5050)

sim_dados_inhomog <- simulate_hawkes_etas_jss(
  T_max = 20,
  X_lim = c(-10, 10), Y_lim = c(-10, 10),
  mu = 0.15, u_fn = u_fn_inhomog, u_max = 1.0,
  A = 0.25, alpha = 2, c_par = 0.2, p_par = 1.3,
  D_par = 0.05, gamma = 1, q_par = 1.5, beta = 2, m_min = 3
)

catalogo_dados_inhomog <- catalog(
    sim_dados_inhomog %>%
    mutate(t = dyears(t),
           data = POSIXct(1) + t,
           long = x, lat = y,
           mag = m,
           date = date(data),
           time = format(data, format = "%H:%M:%S"))
)

plot(catalogo_dados_inhomog)

aj_teste_inhomog_1 <- etas(catalogo_dados_inhomog, nthreads = 8,
                           no.itr = 9)

plot(aj_teste_inhomog_1)

rates(aj_teste_inhomog_1)

# Peru:

set.seed(5050)
sim_dados_homog <- simulate_hawkes_etas_jss(
  T_max = 10,
  X_lim = c(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long)),
  Y_lim = c(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat)),
  mu = 0.005, u_fn = rast_peru, u_max = max(im_peru),
  A = 1.5, alpha = 1.5, c_par = 1.0, p_par = 1.5,
  D_par = 0.2, gamma = 0.1, q_par = 1.1, beta = 2, m_min = 4
)

catalogo_dados_homog <- catalog(
    sim_dados_homog %>%
    mutate(t = dyears(t),
           data = POSIXct(1) + t,
           long = x, lat = y,
           mag = m,
           date = date(data),
           time = format(data, format = "%H:%M:%S"))
)

plot(catalogo_dados_homog)

# Gerando filhos com magnitudes abaixo do valor mínimo:

set.seed(5050)
sim_dados_homog_1 <- simulate_hawkes_etas_zhuang(
  T_max = 10,
  X_lim = c(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long)),
  Y_lim = c(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat)),
  mu = 0.01, u_fn = rast_peru, u_max = max(im_peru),
  A = 0.5, alpha = 1, c_par = 0.2, p_par = 1.15,
  d_par = 0.002, beta = 2, m_min = 4
)

catalogo_dados_homog_1 <- catalog(
    sim_dados_homog_1 %>%
    mutate(t = dyears(t),
           data = POSIXct(1) + t,
           long = x, lat = y,
           mag = m,
           date = date(data),
           time = format(data, format = "%H:%M:%S"))
)

plot(catalogo_dados_homog_1)

#----Gráficos----

ggplot(sim_dados_homog, aes(x = x, y = y)) +
    geom_point(aes(color = gen, size = m), alpha = 0.7) +
    scale_color_viridis_c(option = "plasma", name = "Tempo") +
    scale_size_continuous(name = "Magnitude") +
    theme_minimal() +
    labs(title = "Simulação Hawkes Espaço-Temporal Marcado",
         subtitle = "Kernel baseado em Ogata 1998") +
    coord_fixed()

ggplot(sim_dados_inhomog %>% tail(1000), aes(x = x, y = y)) +
    geom_point(aes(color = gen, size = m), alpha = 0.7) +
    scale_color_viridis_c(option = "plasma", name = "Tempo") +
    scale_size_continuous(name = "Magnitude") +
    theme_minimal() +
    labs(title = "Simulação Hawkes Espaço-Temporal Marcado",
         subtitle = "Kernel baseado em Ogata 1998") +
    coord_fixed()

# Definir os parâmetros do modelo em vetor nomeado:
theta <- c(A = 0.7, alpha = 0.1, c = 0.1, p = 1.15, D = 0.05, gamma = 0.05, q = 2.1, beta = 4)

# Simulação usando função simetas:
set.seed(123)
catalogo_simulado <- simetas(theta, lat.range = c(-5, 5), long.range = c(-5, 5),
                             sim.start = 0, sim.length = 1000, mag.threshold = 4,
                             bkgd = list(x = seq(-5, 5, length.out = 100),
                                         y = seq(-5, 5, length.out = 100),
                                         bkgd = matrix(apply(expand.grid(seq(-5, 5, length.out = 100),
                                                                         seq(-5, 5, length.out = 100)),
                                                             # 1, function(x){(0.01) * u_fn_homog(x[1], x[2])}),
                                                             1, function(x){(0.01) * u_fn_inhomog(x[1], x[2])}),
                                                       nrow = 100)))

# Visualizar o resultado
summary(catalogo_simulado)

plot(catalogo_simulado)

aj_teste_homog_2 <- etas(catalogo_dados_inhomog, nthreads = 8,
                         param0 = c(A = 0.5, alpha = 1.5, c = 0.1, p = 1.15, D = 0.1, gamma = 0.02, q = 2.1, beta = 4))

plot(aj_teste_inhomog_2)

rates(aj_teste_inhomog_2)

#----Projeção 3D----

# Coordenadas x e y do Processo
x <- catalog(iran.quakes)$X$data$x
y <- catalog(iran.quakes)$X$data$y
x <- catalogo_simulado$X$data$x
y <- catalogo_simulado$X$data$y

# Projeção Estereográfica Inversa:
r <- 10
d <- r^2 + x^2 + y^2
X <- 2 * r^2 * x / d
Y <- 2 * r^2 * y / d
Z <- r * ((x^2 + y^2 - r^2) / d)

# Gráfico:
scatterplot3d::scatterplot3d(X, Y, Z, pch = 16,
              highlight.3d = TRUE, main = "Hawkes Projetado na Esfera",
              angle = 30, xlim = c(-r, r), ylim = c(-r, r), zlim = c(-r, r))

library(rgl)

# Interativo:
close3d()
open3d()
plot3d(X, Y, Z,
       xlim = c(-r, r), ylim = c(-r, r), zlim = c(-r, r),
       col = "blue")
points3d(x = x, y = y, z = -r, col = "red")
rglwidget()

#----Simulação Perú----

r <- terra::rast(nrows = 128, ncols = 128,
                 xmin = min(peru.quakes$longlat.coord$long),
                 xmax = max(peru.quakes$longlat.coord$long),
                 ymin = min(peru.quakes$longlat.coord$lat),
                 ymax = max(peru.quakes$longlat.coord$lat))

values(r) <- 1

im_peru <- ppp(x = peru$long, y = peru$lat, window = owin(yrange = c(-23.616257, -1.266735),
                                               xrange = c(-83.05279, -66.82162)),
    marks = peru$mag) %>%
    density.ppp(sigma = bw.scott)

plot(im_peru)
contour(im_peru, add = T)

lat_peru <- range(im_peru$yrow)
long_peru <- range(im_peru$xcol)

theta_peru <- c(A = 0.1, alpha = 0.1, c = 0.1, p = 1.5, D = 0.01, gamma = 0.01, q = 1.001, beta = 4)

set.seed(123)
simu_peru <- simetas(theta_peru, lat.range = lat_peru, long.range = long_peru,
        sim.start = 0, sim.length = 2000, mag.threshold = 4,
        bkgd = list(x = im_peru$xcol,
                    y = im_peru$yrow,
                    bkgd = 0.02 * round(t(im_peru$v))/max(im_peru$v))
        )

plot(simu_peru)

plot(peru.quakes)

debugonce(etas)

etas(simu_peru, nthreads = 12)

param0 <- c(mu = nrow(simu_peru$revents)/(4 * diff(simu_peru$rtperiod) * spatstat.geom::area.owin(simu_peru$region.win)),
            A = 0.5, c = 0.5, alpha = 1, p = 1.1,
            D = 0.5, q = 1.5, gamma = 0.1)

revents <- simu_peru$revents
for (i in 1:1) {
    bkgd <- decluster(param0, spatstat.geom::nndist.default(revents[, 2],
                                                            revents[, 3], k = 5),
                      revents, simu_peru$rpoly, simu_peru$rtperiod, 1000, T, 1)
    revents <- bkgd$revents
}

plot(revents[, 1], revents[, 7], xlab = "time", ylab = "probability of being a background event")

grid <- expand.grid(A = seq(0.1, 1, length.out = 10),
                    c = seq(0.1, 1, length.out = 10))

veros <- numeric(nrow(grid))

for (i in 1:nrow(grid)) {
    print(i)
    theta <- c(mu = 0.002, A = grid[i, 1], c = grid[i, 2],
               alpha = 1, p = 1.1,
               D = 0.5, q = 1.5, gamma = 0.1)
    resultado <- try({
        etasfit(theta = theta, revents = bkgd$revents, rpoly = simu_peru$rpoly,
                tperiod = simu_peru$rtperiod, integ0 = bkgd$integ0,
                ihess = diag(1, 8, 8), verbose = TRUE, ndiv = 1000L,
                eps = 1e-06, cxxcode = TRUE, nthreads = 8L, mver = 1L)
    }, silent = F)
    if (inherits(resultado, "try-error")) {
        veros[i] <- NA
    } else {
        veros[i] <- resultado$loglik
    }
}

ggplot(grid %>% mutate(veros = veros,
                       veros = case_when(near(veros, 0) ~ NA,
                                         .default = veros)),
       aes(x = A, y = c, fill = veros)) +
    geom_tile()

which(veros == max(veros, na.rm = T))
grid[23, ]

ggplot(grid %>% mutate(veros = veros,
                       veros = case_when(near(veros, 0) ~ NA,
                                         veros > 0 ~ NA,
                                         .default = veros)),
       aes(x = A, y = c, fill = veros)) +
    geom_tile()
