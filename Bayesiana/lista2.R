#----Comparação do EQM----

# Gráfico para comparar o Erro Quadrático Médio entre os estimadores:

func_grafico_eqm <- function(n, a, b) {
    m <- cbind(a, b)
    if (nrow(m) == 1) {
        theta <- seq(0, 1, length.out = 1000)
        eqm_mle <- (theta * (1 - theta))/n
        eqm_eap <- (n * theta * (1 - theta) + (a - theta * (a + b))^2)/((n + a + b)^2)
        eqm_map <- ((n * theta * (1 - theta)) + (a - 1 - theta * (a + b - 2))^(2))/((n + a + b - 2)^(2))
        plot(theta, eqm_mle, type = "l", xlab = expression(theta), ylab = "EQM")
        lines(theta, eqm_eap, type = "l", lty = 2, col = "red", lwd = 2)
        lines(theta, eqm_map, type = "l", lty = 3, col = "blue", lwd = 2)
    }
    else {
        par(mfrow = c(round(nrow(m)/2), 2))
        theta <- seq(0, 1, length.out = 1000)
        eqm_mle <- (theta * (1 - theta))/n
        for (i in 1:nrow(m)) {
            eqm_eap <- (n * theta * (1 - theta) + (m[i, 1] - theta * (m[i, 1] + m[i, 2]))^2)/((n + m[i, 1] + m[i, 2])^2)
            eqm_map <- ((n * theta * (1 - theta)) + (m[i, 1] - 1 - theta * (m[i, 1] + m[i, 2] - 2))^(2))/((n + m[i, 1] + m[i, 2] - 2)^(2))
            plot(theta, eqm_mle, type = "l", xlab = expression(theta), ylab = "EQM", main = bquote(a == .(m[i, 1]) ~~~ b == .(m[i, 2])))
            lines(theta, eqm_eap, type = "l", lty = 2, col = "red", lwd = 2)
            lines(theta, eqm_map, type = "l", lty = 3, col = "blue", lwd = 2)
        }
    }
}

func_grafico_eqm(100, rep(c(10, 2), 3), rep(4:6, each = 2))

par(mfrow = c(1, 1))

#----Gráficos de superfície----

# Funções auxiliares:

# EQM do estimador de máxima verossimilhança:
eqm_emv<-function(theta, n){
  theta * (1 - theta)/n
}

# EQM da esperança a posteriori, considerando uma priori Beta(a, b):
eqm_eap <- function(theta, n, a, b){
  (n * theta * (1 - theta) + (a - theta * (a + b))^2)/((n + a + b)^2)
}

# EQM da moda a posteriori, considerando uma priori Beta(a, b):
eqm_map <- function(theta, n, a, b){
    ((n * theta * (1 - theta)) + (a - 1 - theta * (a + b - 2))^(2))/((n + a + b - 2)^(2))
}

# Função para criar o gráfico de superfície:
func_superficie <- function(alpha, beta, n_range, estimador, z_lim = NULL) {
    a <- alpha
    b <- beta
    theta <- seq(0.001, 0.999, length.out = 100)
    n <- seq(n_range[1], n_range[2], length.out = 100)
    z <- switch(estimador,
                "emv" = {outer(theta, n, eqm_emv)},
                "eap" = {outer(theta, n, function(theta, n)(eqm_eap(theta, n, a, b)))},
                "map" = {outer(theta, n, function(theta, n)(eqm_map(theta, n, a, b)))})
    if(is.null(z_lim)) {
        persp(theta, n , z, theta = 45, phi = 20, expand = 1, ticktype = "detailed", shade = 0.75,
              col = 4, zlab = "EQM", main = paste(ifelse(estimador == "emv", "Estimador de máxima verossimilhança",
                                                  ifelse(estimador == "eap", "Esperança a posteriori",
                                                         "Moda a posteriori")), ", a = ", a, ", b = ", b))
    }
    else {
        persp(theta, n , z, theta = 45, phi = 20, expand = 1, ticktype = "detailed", shade = 0.75,
              col = 4, zlab = "EQM", main = paste(ifelse(estimador == "emv", "Estimador de máxima verossimilhança",
                                                  ifelse(estimador == "eap", "Esperança a posteriori",
                                                         "Moda a posteriori")), ", a = ", a, ", b = ", b), zlim = z_lim)
    }
}

# Exemplos de uso:

## Esperança a posteriori, priori = Beta(2, 1):
func_superficie(alpha = 2, beta = 1, n_range = c(1, 10), "eap")

## Comparação do EQM dos três estimadores:
par(mfrow = c(1, 3))
func_superficie(alpha = 2, beta = 1, n_range = c(1, 10), "emv", c(0, 0.25))
func_superficie(alpha = 2, beta = 1, n_range = c(1, 10), "eap", c(0, 0.25))
func_superficie(alpha = 2, beta = 1, n_range = c(1, 10), "map", c(0, 0.25))


#----Outros métodos----

# Para comparar todos em um único gráfico, utilizaremos o pacote plot3D:
library(plot3D)

# Função que plota as três superfícies em um mesmo gráfico:
func_superficie_mult <- function(alpha, beta, n_range) {
    theta <- seq(0.001, 0.999, length.out = 200)
    M <- mesh(theta, seq(n_range[1], n_range[2], length.out = 200))
    z_emv <- outer(theta, seq(n_range[1], n_range[2], length.out = 200), eqm_emv)
    z_eap <- outer(theta, seq(n_range[1], n_range[2], length.out = 200), function(theta, n)(eqm_eap(theta, n, alpha, beta)))
    z_map <- outer(theta, seq(n_range[1], n_range[2], length.out = 200), function(theta, n)(eqm_map(theta, n, alpha, beta)))
    surf3D(M$x, M$y, z_emv, clim = c(0, 0.25), phi = 20, theta = 45, border = "black", bty = "f")
    surf3D(M$x, M$y, z_eap, clim = c(0, 0.25), add = T, border = "gray")
    surf3D(M$x, M$y, z_map, clim = c(0, 0.25), add = T)
}

# Exemplo de uso:
func_superficie_mult(2, 4, c(1, 10))

# A superfície que está com border preto é a de verossimilhança,
# a em cinza claro é a esperança a posteriori e a sem border é a moda a posteriori

# Também é possível tornar a superfície interativa, com o pacote plotly:
library(plotly)

# Função que plota as três superfícies em um mesmo gráfico interativo:
func_superficie_mult_plotly <- function(alpha, beta, n_range) {
    theta <- seq(0.001, 0.999, length.out = 200)
    M <- mesh(theta, seq(n_range[1], n_range[2], length.out = 200))
    z_emv <- outer(theta, seq(n_range[1], n_range[2], length.out = 200), eqm_emv)
    z_eap <- outer(theta, seq(n_range[1], n_range[2], length.out = 200), function(theta, n)(eqm_eap(theta, n, alpha, beta)))
    z_map <- outer(theta, seq(n_range[1], n_range[2], length.out = 200), function(theta, n)(eqm_map(theta, n, alpha, beta)))
    plot_ly(x = theta, y = seq(n_range[1], n_range[2], length.out = 200), cmin = 0, cmax = 0.25) %>%
        add_surface(z = z_emv, contours = list(x = list(show = T, size = 0.01, color = "black"),
                                               y = list(show = T, size = 0.01, color = "black"),
                                               z = list(show = T, size = 0.01, color = "black"))) %>%
        add_surface(z = z_eap, contours = list(x = list(show = T, size = 0.01, color = "white"),
                                               y = list(show = T, size = 0.01, color = "white"),
                                               z = list(show = T, size = 0.01, color = "white"))) %>%
        add_surface(z = z_map)
}

# Exemplo de uso:
func_superficie_mult_plotly(2, 4, c(1, 20))

# Apesar de serem interativos, os gráficos feitos com plot_ly são mais
# difíceis de visualizar, por causa das paletas de cores, que não consegui alterar.
