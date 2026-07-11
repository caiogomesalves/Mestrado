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
  T_max = 1000,
  X_lim = c(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long)),
  Y_lim = c(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat)),
  mu = mu_global,
  u_fn = terra::rast(rasters[[length(rasters)]]%>%as.im(W=owin(xrange = c(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long)),
                                                               yrange = c(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat))))),
  u_max = max(rasters[[length(rasters)]]),
  A = thetas[[length(thetas)]]["A"],
  alpha = thetas[[length(thetas)]]["alpha"],
  c_par = thetas[[length(thetas)]]["c"],
  p_par = thetas[[length(thetas)]]["p"],
  d_par = thetas[[length(thetas)]]["D"],
  beta = beta_GR, m_min = 4
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

#----Teste----

theta_estimado <- thetas[[length(thetas)]]

set.seed(5050)
dados_simulados <- simulate_hawkes_etas_zhuang(
    T_max = max(peru.quakes$rtperiod),
    X_lim = range(peru.quakes$region.poly$long),
    Y_lim = range(peru.quakes$region.poly$lat),
    mu = mu_global,
    u_fn = terra::rast(rasters[[length(rasters)]] %>%
                       as.im(W = owin(xrange = range(peru.quakes$region.poly$long),
                                      yrange = range(peru.quakes$region.poly$lat)))),
    u_max = max(rasters[[length(rasters)]]),
    A = theta_estimado["A"], alpha = theta_estimado["alpha"],
    c_par = theta_estimado["c"], p_par = theta_estimado["p"],
    d_par = theta_estimado["D"],
    beta = 2, m_min = 4
)

catalogo_dados_simulados <- catalog(
    dados_simulados %>%
    mutate(t = ddays(t),
           data = POSIXct(1) + t,
           long = x, lat = y,
           mag = m,
           date = date(data),
           time = format(data, format = "%H:%M:%S"))
)

plot(catalogo_dados_simulados)

lambda_total_bg <- sum(rasters[[length(rasters)]])

# Parâmetros e sementes prévias...
beta_GR <- nrow(peru)/sum((peru$mag - 4))

# 1. Simular Background
n_bg <- rpois(1, lambda_total_bg)
df_sim <- data.frame(
  id = 1:n_bg,
  parent_id = 0, # 0 indica background
  generation = 0,
  time = runif(n_bg, 0, T_max),
  # sample_locations() é uma função sua para sortear (x,y) pelo raster de probabilidades
  x = ...,
  y = ...,
  mag = M0 + rexp(n_bg, rate = beta_GR)
)

# Fila de pais para a próxima iteração
pais_atuais <- df_sim

# 2. Simular Ramificações
geracao_atual <- 1
id_counter <- n_bg

while(nrow(pais_atuais) > 0) {

  # Quantos filhos CADA pai produzirá?
  n_filhos <- rpois(nrow(pais_atuais), A * exp(alpha * (pais_atuais$mag - M0)))

  # Filtra apenas os pais que tiveram ao menos 1 filho
  pais_férteis <- pais_atuais[n_filhos > 0, ]
  n_filhos <- n_filhos[n_filhos > 0]

  if(length(n_filhos) == 0) break # Extinção atingida

  # Expande as linhas dos pais férteis para parear com cada filho
  # ex: se o pai 3 teve 2 filhos, a linha dele é duplicada
  df_filhos <- pais_férteis[rep(seq_len(nrow(pais_férteis)), n_filhos), ]

  # Gera deltas temporais (Omori inversa)
  u <- runif(nrow(df_filhos))
  dt <- c * (u^(-1 / (p - 1)) - 1)
  df_filhos$time <- df_filhos$time + dt

  # Gera deltas espaciais (Normal dependente da magnitude)
  sigma <- sqrt(D * exp(alpha * (df_filhos$mag - M0)))
  dx <- rnorm(nrow(df_filhos), mean = 0, sd = sigma)
  dy <- rnorm(nrow(df_filhos), mean = 0, sd = sigma)
  df_filhos$x <- df_filhos$x + dx
  df_filhos$y <- df_filhos$y + dy

  # Novas magnitudes
  df_filhos$mag <- M0 + rexp(nrow(df_filhos), rate = beta_GR)

  # Atualiza metadados
  df_filhos$parent_id <- df_filhos$id
  df_filhos$id <- (id_counter + 1):(id_counter + nrow(df_filhos))
  df_filhos$generation <- geracao_atual

  # FILTRO: Mantém apenas filhos dentro da caixa de tempo e espaço
  df_filhos <- subset(df_filhos, time <= T_max & x >= x_min & x <= x_max & y >= y_min & y <= y_max)

  # Acumula no catálogo principal
  df_sim <- rbind(df_sim, df_filhos)

  # Define os filhos sobreviventes como a próxima geração de pais
  pais_atuais <- df_filhos
  id_counter <- max(df_sim$id)
  geracao_atual <- geracao_atual + 1
}

# Ordena o catálogo final por tempo de ocorrência
df_sim <- df_sim[order(df_sim$time), ]

#====Simulador Pareto====
# =============================================================================
# Simulação de catálogo ETAS espaço-temporal (kernel de Pareto)
# via representação de ramificação (cluster representation de Hawkes & Oakes, 1974)
#
# Ideia central:
#   Cada evento i (background ou réplica) gera um número de "filhos" diretos
#   que é Poisson(kappa(m_i)), pois g(.) e f(.) sao densidades normalizadas
#   (integram a 1 sobre [0,Inf) e R^2, respectivamente). Isso e' exato --
#   nao e' uma aproximacao por thinning generico, e' a representacao teorica
#   do proprio processo (Hawkes & Oakes, 1974; Zhuang, Ogata & Vere-Jones, 2004).
# =============================================================================


# -----------------------------------------------------------------------------
# Checagem de criticidade (razao de ramificacao esperada)
# n_bar = A * beta / (beta - alpha), com beta = b*log(10) (Gutenberg-Richter)
# Precisa de beta > alpha para a razao ser finita, e n_bar < 1 para o processo
# ser subcritico (catalogo nao explode, media de eventos finita).
# Ver Helmstetter & Sornette (2002) para essa formula da razao de ramificacao.
# -----------------------------------------------------------------------------
verificar_criticidade <- function(theta, b) {
  A <- unname(theta["A"]); alpha <- unname(theta["alpha"])
  beta <- b * log(10)
  if (beta <= alpha) {
    warning(sprintf(
      "beta (%.3f) <= alpha (%.3f): razao de ramificacao e' infinita/mal definida. ",
      beta, alpha),
      "A distribuicao de magnitudes tem cauda mais pesada que o decaimento ",
      "da produtividade -- catalogo pode explodir numericamente.")
    return(invisible(NULL))
  }
  n_bar <- A * beta / (beta - alpha)
  msg <- sprintf("Razao de ramificacao estimada: n = %.4f", n_bar)
  if (n_bar >= 1) {
    warning(msg, " -- processo SUPERCRITICO (n >= 1). ",
            "O numero esperado de eventos e' infinito/instavel; a simulacao ",
            "provavelmente vai bater no limite de geracoes ou demorar muito.")
  } else {
    message(msg, " (subcritico, processo estavel)")
  }
  invisible(n_bar)
}


# -----------------------------------------------------------------------------
# Funcao principal de simulacao
# -----------------------------------------------------------------------------
simular_catalogo_etas <- function(
                                  theta,                          # vetor nomeado: A, alpha, c, p, D, q, gamma
                                  M0      = 4.5,                  # magnitude de corte (completude)
                                  b       = 1.0,                  # valor-b de Gutenberg-Richter p/ magnitudes
                                  Mmax    = NULL,                 # Magnitude maxima (NULL = sem truncamento)
                                  T_max,                          # Janela temporal de simulacao [0, T_max]
                                  mu_type = c("constant", "raster"),
                                  mu_const  = NULL,               # Usado se mu_type == "constant" (eventos/dia/grau^2)
                                  mu_raster = NULL,               # Usado se mu_type == "raster" (matriz nx x ny)
                                  x_grid = NULL, y_grid = NULL,   # Necessarios se mu_type == "raster"
                                  xlim = NULL, ylim = NULL,       # Janela espacial (usada se mu_type == "constant")
                                  max_geracoes = 75,              # Trava de seguranca contra explosão
                                  checar_criticidade = TRUE,
                                  semente = NULL
                                  ) {
    if (!is.null(semente)){
        set.seed(semente)
    }
    mu_type <- match.arg(mu_type)
    if (checar_criticidade){
        verificar_criticidade(theta, b)
    }
    A     <- unname(theta["A"])
    alpha <- unname(theta["alpha"])
    cc    <- unname(theta["c"])
    pp    <- unname(theta["p"])
    D     <- unname(theta["D"])
    q     <- unname(theta["q"])
    gamma <- unname(theta["gamma"])
    # Garante que valores sejam próprios para densidades:
    stopifnot(pp > 1, q > 1)
    # ---- Magnitudes: Gutenberg-Richter (truncada se Mmax fornecido) ----
    simular_magnitudes <- function(n) {
        if (n == 0) return(numeric(0))
        beta <- b * log(10)
        if (is.null(Mmax)) {
            M0 - log(runif(n)) / beta
        } else {
            u <- runif(n)
            Fmax <- 1 - exp(-beta * (Mmax - M0))
            M0 - log(1 - u * Fmax) / beta
        }
    }
    # ---- Localizacao de background a partir de um raster (grid) ----
    amostrar_local_raster <- function(n) {
        # Caso nenhuma observação retorna data frame vazio:
        if (n == 0){
            return(data.frame(x = numeric(0), y = numeric(0)))
        }
        n_x <- length(x_grid)
        n_y <- length(y_grid)
        probs <- as.vector(mu_raster)
        probs <- probs / sum(probs)
        idx <- sample.int(length(probs), size = n, replace = TRUE, prob = probs)
        col_idx <- ((idx - 1) %/% n_x) + 1
        row_idx <- ((idx - 1) %% n_x) + 1
        dx <- diff(x_grid)[1]; dy <- diff(y_grid)[1]
        x <- x_grid[row_idx] + runif(n, -dx / 2, dx / 2)
        y <- y_grid[col_idx] + runif(n, -dy / 2, dy / 2)
        data.frame(x = x, y = y)
    }
    # ---- Deslocamento espacial de replicas: kernel de Pareto isotropico ----
    # CDF do raio (obtida integrando f(x,y) em coordenadas polares):
    #   F(r) = 1 - (1 + r^2/S)^{-(q-1)}   =>   r^2 = S * [(1-u)^{-1/(q-1)} - 1]
    simular_offset_espacial <- function(mag) {
        n <- length(mag)
        if (n == 0) {
            return(data.frame(dx = numeric(0), dy = numeric(0)))
        }
        S_i <- D * exp(gamma * (mag - M0))
        u   <- runif(n)
        r2  <- S_i * ((1 - u)^(-1 / (q - 1)) - 1)
        r   <- sqrt(r2)
        ang <- runif(n, 0, 2 * pi)
        data.frame(dx = r * cos(ang), dy = r * sin(ang))
    }
    # ---- Deslocamento temporal: Omori-Utsu (CDF inversa) ----
    # F(t) = 1 - (c/(t+c))^{p-1}   =>   t = c * [(1-u)^{-1/(p-1)} - 1]
    simular_offset_temporal <- function(n) {
        if (n == 0) return(numeric(0))
        u <- runif(n)
        cc * ((1 - u)^(-1 / (pp - 1)) - 1)
    }
    # ---- 1. Eventos de background (geracao 0) ----
    if (mu_type == "constant") {
        stopifnot(!is.null(mu_const), !is.null(xlim), !is.null(ylim))
        area <- diff(xlim) * diff(ylim)
        n_bg <- rpois(1, mu_const * area * T_max)
        x0 <- runif(n_bg, xlim[1], xlim[2])
        y0 <- runif(n_bg, ylim[1], ylim[2])
    } else {
        stopifnot(!is.null(mu_raster), !is.null(x_grid), !is.null(y_grid))
        dxg <- diff(x_grid)[1]
        dyg <- diff(y_grid)[1]
        # Integral aproximada de mu(x,y):
        mu_total <- sum(mu_raster) * dxg * dyg
        n_bg <- rpois(1, mu_total * T_max)
        loc <- amostrar_local_raster(n_bg)
        x0 <- loc$x
        y0 <- loc$y
    }
    t0 <- runif(n_bg, 0, T_max)
    m0 <- simular_magnitudes(n_bg)
    catalogo <- data.frame(
        id = seq_len(n_bg), parent = 0L, geracao = 0L,
        time = t0, x = x0, y = y0, mag = m0
    )
    # ---- 2. Processar réplicas geração por geração (fila / BFS) ----
    proximo_id <- n_bg + 1L
    fila <- catalogo
    g_final <- 0L
    if (n_bg > 0) {
        for (g in seq_len(max_geracoes)) {
            g_final <- g
            if (nrow(fila) == 0){
                break
            }
            kappa_vals <- A * exp(alpha * (fila$mag - M0))
            n_filhos   <- rpois(nrow(fila), kappa_vals)
            if (sum(n_filhos) == 0){
                break
            }
            idx_pais <- rep(seq_len(nrow(fila)), times = n_filhos)
            pais <- fila[idx_pais, , drop = FALSE]
            mags_filhos <- simular_magnitudes(nrow(pais))
            dt          <- simular_offset_temporal(nrow(pais))
            off_esp     <- simular_offset_espacial(pais$mag)
            filhos <- data.frame(
                id      = proximo_id + seq_len(nrow(pais)) - 1L,
                parent  = pais$id,
                geracao = g,
                time    = pais$time + dt,
                x       = pais$x + off_esp$dx,
                y       = pais$y + off_esp$dy,
                mag     = mags_filhos
            )
            proximo_id <- proximo_id + nrow(filhos)
            # Descarta filhos fora da janela temporal
            filhos <- filhos[filhos$time <= T_max, , drop = FALSE]
            catalogo <- rbind(catalogo, filhos)
            fila <- filhos
        }
        if (g_final == max_geracoes && nrow(fila) > 0) {
            warning("Número máximo de gerações atingido -- verifique a razão de ",
                    "ramificação (verificar_criticidade) antes de confiar no catálogo.")
        }
    }
    catalogo <- catalogo[order(catalogo$time), ]
    rownames(catalogo) <- NULL
    attr(catalogo, "theta_verdadeiro") <- theta
    attr(catalogo, "M0") <- M0
    attr(catalogo, "b")  <- b
    catalogo
}

# =============================================================================
# Exemplo de uso: estudo de recuperacao de parametros
# =============================================================================
if (FALSE) {

  theta_verdadeiro <- c(
    A = 0.05, alpha = 1.8,
    c = 0.01, p = 1.15,
    D = 0.02, q = 1.6, gamma = 1.0
  )

  verificar_criticidade(theta_verdadeiro, b = 1.0)

  set.seed(123)
  cat_sim <- simular_catalogo_etas(
    theta   = thetas[[length(thetas)]],
    M0      = 4.5,
    b       = 1.0,
    Mmax    = 10,
    T_max   = 20000,                 # 10 anos, por exemplo
    mu_type = "raster",
    mu_raster = rasters[[length(rasters)]],
    x_grid = x_grid, y_grid = y_grid,
    # mu_const = 0.002,                # ajuste para gerar um numero razoavel de eventos
    # xlim = c(-80, -68), ylim = c(-18, 0),   # aprox. janela do Peru
    semente = 5050
  )

  nrow(cat_sim)
  table(cat_sim$geracao)            # quantos eventos por geracao (deve decair)
  mean(cat_sim$parent == 0)         # fracao de background observada

  # Depois: rode SUA rotina de ajuste (thinning_pareto_zurich.R) em cat_sim
  # e compare theta_estimado vs theta_verdadeiro -- idealmente repetindo
  # varias vezes (ex: 100-500 replicas) para ver vies e variancia do estimador,
  # nao so' uma unica simulacao.

  # Para usar o SEU raster de background real (mais realista, reproduz
  # inclusive os efeitos de borda que discutimos):
  #
  # cat_sim2 <- simular_catalogo_etas(
  #   theta = theta_verdadeiro, M0 = 4.5, b = 1.0, T_max = 3650,
  #   mu_type = "raster", mu_raster = rasters[[length(rasters)]],
  #   x_grid = x_grid, y_grid = y_grid
  # )
}

#====Simulações Irã====
thetas_ajustados_iran <- list()
backgrounds_iran <- list()

for (k in 1:10) {
    cat_sim_iran <- simular_catalogo_etas(
        theta   = thetas_iran[[length(thetas_iran)]],
        M0      = 4.5,
        b       = 2.0,
        Mmax    = 10,
        T_max   = 3650,
        mu_type = "raster",
        mu_raster = rasters_iran[[length(rasters_iran)]],
        x_grid = x_grid_iran,
        y_grid = y_grid_iran,
        semente = k
    )
    cat_temp_iran <- catalog(
        cat_sim_iran %>%
        mutate(t = ddays(time),
               data = POSIXct(1) + t,
               long = x, lat = y,
               date = date(data),
               time = format(data, format = "%H:%M:%S")),
        region.poly = iran.cat$region.poly
    )
    theta_iran_simu <- c("A" = 0.23, "alpha" = 2.8, "c" = 0.022,
                         "p" = 1.12, "D" = 0.012, "q" = 2.3, "gamma" = 0.03)
    thetas_iran_temp <- list()
    thetas_iran_temp[[1]] <- theta_iran_simu
    # Número total de eventos:
    n_events_iran_simu <- nrow(cat_temp_iran$revents)
    # Definir as sequências do grid explicitamente
    x_grid_iran <- seq(min(iran.cat$longlat.coord$long), max(iran.cat$longlat.coord$long), length.out = 128)
    y_grid_iran <- seq(min(iran.cat$longlat.coord$lat), max(iran.cat$longlat.coord$lat), length.out = 128)
    # Raster inicial:
    r_iran_temp <- matrix(1, nrow = 128, ncol = 128)
    # Na lista de rasters, guardaremos as matrizes:
    rasters_iran_temp <- list()
    rasters_iran_temp[[1]] <- r_iran_temp
    # Data.frame com os eventos para funções:
    df_temp_iran_temp <- data.frame(x = cat_temp_iran$longlat.coord$long,
                                    y = cat_temp_iran$longlat.coord$lat,
                                    mag = cat_temp_iran$revents[, "mm"] + 4.5,
                                    time = cat_temp_iran$revents[, "tt"])
    df_temp_iran_temp$time <- df_temp_iran_temp$time - min(df_temp_iran_temp$time)
    # Distâncias para kernel adaptativo:
    d_j_iran_temp <- d_j(cat_temp_iran$longlat.coord[, c("long", "lat")])
    num_cores <- 8
    # Paralelização do cálculo de probabilidades:
    probabilidades_iran_temp <- mclapply(1:n_events_iran_simu, function(j) {
        # Primeiro evento tem probabilidade 0 de ser descendente:
        if (j == 1) return(numeric(0))
        # Vetor de índices dos antecessores:
        i_indices <- 1:(j - 1)
        # Denominador comum para cada evento j:
        denom <- lambda_total(
            df_temp_iran_temp$time[j], df_temp_iran_temp$x[j], df_temp_iran_temp$y[j], df_temp_iran_temp$mag[j],
            4.5, df_temp_iran_temp, theta_iran_simu, r_iran_temp, x_grid_iran, y_grid_iran
        )
        # Computação do numerador para todos os antecessores:
        numeradores <- kappa_func(df_temp_iran_temp$mag[i_indices], 4.5, theta_iran_simu["A"], theta_iran_simu["alpha"]) *
            g_func(df_temp_iran_temp$time[j] - df_temp_iran_temp$time[i_indices], theta_iran_simu["c"], theta_iran_simu["p"]) *
            f_spatial(
            (df_temp_iran_temp$x[j] - df_temp_iran_temp$x[i_indices])^2 + (df_temp_iran_temp$y[j] - df_temp_iran_temp$y[i_indices])^2,
            df_temp_iran_temp$mag[i_indices], 4.5, theta_iran_simu["D"], theta_iran_simu["q"], theta_iran_simu["gamma"]
            )
        # Retorna o vetor p_{i,j} para o evento j:
        return(numeradores / denom)
    }, mc.cores = num_cores)
    # Probabilidades totais:
    prob_total_iran_temp <- probabilidades_iran_temp %>%
        lapply(sum) %>%
        do.call(c, .)
    # Estima o background rate:
    r_iran_temp <- estimate_background_kernel(x_grid_iran,
                                              y_grid_iran,
                                              data.frame(x = cat_temp_iran$longlat.coord$long,
                                                         y = cat_temp_iran$longlat.coord$lat,
                                                         rho = prob_total_iran_temp),
                                              d_j_iran_temp,
                                              T_total = max(df_temp_iran_temp$time))
    otimizacao_iran_temp <- optim(theta_iran_simu,
                                  log_lik_etas,
                                  gr = gradiente_etas_completo,
                                  df_eventos = df_temp_iran_temp,
                                  matriz_p_ij = probabilidades_iran_temp,
                                  prob_total = prob_total_iran_temp,
                                  raster_mu = r_iran_temp, x_grid = x_grid_iran, y_grid = y_grid_iran,
                                  method = "L-BFGS-B",
                                  lower = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-6),
                                  upper = c(15, 10, 15, 15, 15, 15, 10),
                                  control = list(maxit = 1000),
                                  hessian = T)
    theta_iran_simu <- otimizacao_iran_temp$par
    thetas_iran_temp[[2]] <- theta_iran_simu
    rasters_iran_temp[[2]] <- r_iran_temp
    i <- 2
    # Considerando tolerância de 10e-6:
    while ((any(abs(rasters_iran_temp[[i]] - rasters_iran_temp[[i - 1]]) > 1e-3)) | (any(abs(thetas_iran_temp[[i]] - thetas_iran_temp[[i - 1]]) > 1e-3))) {
        # Calcula pares de probabilidades:
        probabilidades_iran_temp <- mclapply(1:n_events_iran_simu, function(j) {
            if (j == 1) return(numeric(0))
            i_indices <- 1:(j - 1)
            denom <- lambda_total(
                df_temp_iran_temp$time[j], df_temp_iran_temp$x[j], df_temp_iran_temp$y[j], df_temp_iran_temp$mag[j],
                4.5, df_temp_iran_temp, theta_iran_simu, r_iran_temp, x_grid_iran, y_grid_iran
            )
            numeradores <- kappa_func(df_temp_iran_temp$mag[i_indices], 4.5, theta_iran_simu["A"], theta_iran_simu["alpha"]) *
                g_func(df_temp_iran_temp$time[j] - df_temp_iran_temp$time[i_indices], theta_iran_simu["c"], theta_iran_simu["p"]) *
                f_spatial(
                (df_temp_iran_temp$x[j] - df_temp_iran_temp$x[i_indices])^2 + (df_temp_iran_temp$y[j] - df_temp_iran_temp$y[i_indices])^2,
                df_temp_iran_temp$mag[i_indices], 4.5, theta_iran_simu["D"], theta_iran_simu["q"], theta_iran_simu["gamma"]
                )
            return(numeradores / denom)
        }, mc.cores = num_cores)
        # Calcula probabilidades totais:
        prob_total_iran_temp <- probabilidades_iran_temp %>%
            lapply(sum) %>%
            do.call(c, .)
        # Estima parâmetros theta:
        otimizacao_iran_temp <- optim(theta_iran_simu,
                                      log_lik_etas,
                                      gr = gradiente_etas_completo,
                                      df_eventos = df_temp_iran_temp,
                                      matriz_p_ij = probabilidades_iran_temp,
                                      prob_total = prob_total_iran_temp,
                                      raster_mu = r_iran_temp, x_grid = x_grid_iran, y_grid = y_grid_iran,
                                      method = "L-BFGS-B",
                                      lower = c(1e-6, 1e-6, 1e-6, 1.001, 1e-6, 1.001, 1e-6),
                                      upper = c(15, 10, 15, 15, 15, 15, 10),
                                      control = list(maxit = 1000),
                                      hessian = T)
        theta_iran_simu <- otimizacao_iran_temp$par
        thetas_iran_temp[[i + 1]] <- theta_iran_simu
        # Estima background rate:
        r_iran_temp <- estimate_background_kernel(x_grid_iran,
                                                  y_grid_iran,
                                                  data.frame(x = cat_temp_iran$longlat.coord$long,
                                                             y = cat_temp_iran$longlat.coord$lat,
                                                             rho = prob_total_iran_temp),
                                                  d_j_iran_temp,
                                                  T_total = max(df_temp_iran_temp$time))
        rasters_iran_temp[[i + 1]] <- r_iran_temp
        print(i + 1)
        # Atualiza índice:
        i <- i + 1
    }
    thetas_ajustados_iran[[k]] <- thetas_iran_temp[[length(thetas_iran_temp)]]
    backgrounds_iran[[k]] <- rasters_iran_temp[[length(rasters_iran_temp)]]
}
