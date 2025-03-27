#----Questão 1----

# Seed para reprodutibilidade:
set.seed(5050)

# Valores dos x:
x <- 1:10

# 10.000 simulações dos Y_i:
y_sim <- replicate(10000, 5 + 2 * x + rnorm(10, 0, 1))

# Vetores e listas auxiliares para manter as informações:
beta_0 <- vector(mode = "numeric", 10000)
beta_1 <- vector(mode = "numeric", 10000)
sigma_2 <- vector(mode = "numeric", 10000)

ic_beta_0 <- vector(mode = "list", 10000)
ic_beta_1 <- vector(mode = "list", 10000)

estat_teste_beta_0 <- vector(mode = "numeric", 10000)
estat_teste_beta_1_2 <- vector(mode = "numeric", 10000)
estat_teste_beta_1_1.8 <- vector(mode = "numeric", 10000)

# Função auxiliar:
func_soma <- function(a, b) {
    sum((a - mean(a)) * (b - mean(b)))
}

# Loop para fazer os ajustes e popular os vetores e listas:
for (i in 1:10000) {
    # Estimativas pontuais dos Betas:
    beta_1[i] <- func_soma(y_sim[, i], x)/func_soma(x, x)
    beta_0[i] <- mean(y_sim[, i]) - beta_1[i] * mean(x)
    # Estimativa da variância pelos resíduos:
    sigma_2[i] <- sum((y_sim[, i] - (beta_0[i] +
                                     beta_1[i] * x))^(2))/(length(x) - 2)
    #sigma_2[i] <- sum(residuals(ajuste)^(2))/(length(x) - 2)
    # Valor da distribuição t para criação dos intervalos de confiança:
    t_alpha <- qt(1 - 0.05/2, length(x) - 2)
    # Calcula os intervalos de confiança dos betas:
    ic_beta_0[[i]] <- c(`2.5 %` = beta_0[i] -
                            (sqrt((sigma_2[i] * sum(x^(2)))/
                                  (length(x) * func_soma(x, x))) * t_alpha),
                        `97.5 %` = beta_0[i] +
                            (sqrt((sigma_2[i] * sum(x^(2)))/
                                  (length(x) * func_soma(x, x))) * t_alpha))
    ic_beta_1[[i]] <- c(`2.5 %` = beta_1[i] -
                            sqrt((sigma_2[i])/
                                 (func_soma(x, x))) * t_alpha,
                        `97.5 %` = beta_1[i] +
                            sqrt((sigma_2[i])/
                                 (func_soma(x, x))) * t_alpha)
    # Estatística de teste para Beta_0 = 5:
    estat_teste_beta_0[i] <- (beta_0[i] - 5)/sqrt((sigma_2[i] * sum(x^(2)))/
                                                  (length(x) * func_soma(x, x)))
    # Estatística de teste para Beta_1 = 2:
    estat_teste_beta_1_2[i] <- (beta_1[i] - 2)/sqrt(sigma_2[i]/func_soma(x, x))
    # Estatística de teste para Beta_1 = 1.8:
    estat_teste_beta_1_1.8[i] <- (beta_1[i] - 1.8)/sqrt(sigma_2[i]/func_soma(x, x))
}

# Transformação dos intervalos de confiança em matrizes:
ic_beta_0 <- do.call(rbind, ic_beta_0)
ic_beta_1 <- do.call(rbind, ic_beta_1)

#----Questão 2----

# Pacote para manipular os dados:
library(tidyverse)

data.frame(x = 1:10000, y = beta_0) %>%
    ggplot(aes(x = y)) +
    geom_histogram() +
    geom_vline(xintercept = 5, colour = "red") +
    labs(y = "Contagem", x = expression(hat(beta)[0]))

# Gráfico das estimativas pontuais de Beta_0 ao longo das simulações:
data.frame(x = 1:10000, y = beta_0) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point() +
    geom_hline(yintercept = 5, colour = "red") +
    labs(x = "Simulação", y = expression(hat(beta)[0]))

# Média das estimativas:
mean(beta_0)

# O valor esperado é 5.

# Variância das estimativas:
var(beta_0)

# O valor esperado para a variância é dado por:

var_beta_0_teo <- (1 * sum(x^(2)))/(length(x) * sum((x - mean(x))^(2)))

# Ambos se aproximam muito dos valores esperados.

# Gráfico das estimativas pontuais de Beta_1 ao longo das simulações:
data.frame(x = 1:10000, y = beta_1) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point() +
    geom_hline(yintercept = 2, colour = "red") +
    labs(x = "Simulação", y = expression(hat(beta)[1]))

# Média das estimativas:
mean(beta_1)

# O valor esperado é 2.

# Variância das estimativas:
var(beta_1)

# O valor esperado para a variância é dado por:

var_beta_1_teo <- (1/sum((x - mean(x))^(2)))

# Novamente, os valores se aproximam muito dos valores esperados.

#----Questão 3----

## Para Beta_0:

# Quantidade de intervalos de confiança para Beta_0
# que não contém o valor verdadeiro (Beta_0 = 5):
sum(ic_beta_0[, 1] > 5 | ic_beta_0[, 2] < 5)

# Logo, a quantidade de intervalos que contém é:
10000 - sum(ic_beta_0[, 1] > 5 | ic_beta_0[, 2] < 5)

# Como a função confint, por padrão, constrói intervalos de
# 95% de confiança, esperava-se que 9500 simulações gerassem
# intervalos que compreendam o valor verdadeiro para Beta_0

# Podemos visualizar os 200 primeiros intervalos
# a fim de verificar quantos contém o verdadeiro valor de Beta_0:
as.data.frame(ic_beta_0) %>%
    rownames_to_column() %>%
    mutate(rowname = as.numeric(rowname),
           Intervalo = ifelse(`2.5 %` > 5 | `97.5 %` < 5, "Não Contém", "Contém")) %>%
    filter(rowname <= 200) %>%
    ggplot(aes(x = rowname, y = `2.5 %`)) +
    geom_segment(aes(yend = `97.5 %`, colour = Intervalo)) +
    geom_hline(yintercept = 5, colour = "blue") +
    labs(x = "Simulação", y = expression(hat(beta)[0])) +
    scale_colour_manual(values = c("black", "red")) +
    theme_bw()

# 11 dos 200 intervalos não contém o valor Beta_0 = 5
# (que são os marcados em vermelho), que é aproximadamente 5%

## Para Beta_1:

# Quantidade de intervalos de confiança para Beta_1
# que não contém o valor verdadeiro (Beta_1 = 2):
sum(ic_beta_1[, 1] > 2 | ic_beta_1[, 2] < 2)

# Logo, a quantidade de intervalos que contém é:
10000 - sum(ic_beta_1[, 1] > 2 | ic_beta_1[, 2] < 2)

# Como a função confint, por padrão, constrói intervalos de
# 95% de confiança, esperava-se que 9500 simulações gerassem
# intervalos que compreendam o valor verdadeiro para Beta_1

# Podemos visualizar os 200 primeiros intervalos
# a fim de verificar quantos contém o verdadeiro valor de Beta_1:
as.data.frame(ic_beta_1) %>%
    rownames_to_column() %>%
    mutate(rowname = as.numeric(rowname),
           Intervalo = ifelse(`2.5 %` > 2 | `97.5 %` < 2, "Não Contém", "Contém")) %>%
    filter(rowname <= 200) %>%
    ggplot(aes(x = rowname, y = `2.5 %`)) +
    geom_segment(aes(yend = `97.5 %`, colour = Intervalo)) +
    geom_hline(yintercept = 2, colour = "blue") +
    labs(x = "Simulação", y = expression(hat(beta)[0])) +
    scale_colour_manual(values = c("black", "red")) +
    theme_bw()

# 13 dos 200 intervalos não contém o valor Beta_1 = 2
# (que são os marcados em vermelho), que é aproximadamente 5%

## Para ambos os Betas:

# Quantidade de simulações em que ambos os intervalos
# não contém os valores verdadeiros para Beta_0 e Beta_1:
(sum((ic_beta_0[, 1] > 5 | ic_beta_0[, 2] < 5) &
    (ic_beta_1[, 1] > 2 | ic_beta_1[, 2] < 2)))

# Logo, a quantidade de intervalos que contém é:
10000 - sum((ic_beta_0[, 1] > 5 | ic_beta_0[, 2] < 5) &
            (ic_beta_1[, 1] > 2 | ic_beta_1[, 2] < 2))

# TODO: explicar a covariância entre os Betas e a quantidade de simulações
# em que ambos os intervalos contém os valores verdadeiros dos Betas.

#----Questão 4----

# Comparação do valor absoluto das estatísticas de teste obtidas com
# o quantil (1-0.05/2)% da distribuição t com 8 graus de liberdade:
sum(abs(estat_teste_beta_0) > qt(1 - 0.05/2, length(x) - 2))

# Aqui temos 498 simulações em que a hipótese seria rejeitada,
# que condiz com a quantidade de intervalos de confiança que não
# contém o verdadeiro valor de Beta_0.

# Esperava-se que a quantidade de rejeições fosse de 5%, ou seja,
# 500 das simulações.

#----Questão 5----

# Comparação do valor absoluto das estatísticas de teste obtidas com
# o quantil (1-0.05/2)% da distribuição t com 8 graus de liberdade:

sum(abs(estat_teste_beta_1_2) > qt(1 - 0.05/2, length(x) - 2))

# Aqui temos 534 simulações em que a hipótese seria rejeitada,
# que condiz com a quantidade de intervalos de confiança que não
# contém o verdadeiro valor de Beta_1.

# Esperava-se que a quantidade de rejeições fosse de 5%, ou seja,
# 500 das simulações.

#----Questão 6----

# Comparação do valor absoluto das estatísticas de teste obtidas com
# o quantil (1-0.05/2)% da distribuição t com 8 graus de liberdade:

sum(abs(estat_teste_beta_1_1.8) > qt(1 - 0.05/2, length(x) - 2))

# Aqui temos 534 simulações em que a hipótese seria rejeitada,
# que condiz com a quantidade de intervalos de confiança que não
# contém o valor de Beta_1 = 1.8:

sum(ic_beta_1[, 1] > 1.8 | ic_beta_1[, 2] < 1.8)

# A quantidade esperada é correspondente ao valor da diferença entre
# o valor médio de Beta_1 (que é aproximadamente 2) para 1.8 na distribuição
# t com 8 graus de liberdade vezes a quantidade de simulações, que é dado por:
dt(1.8 - mean(beta_1), length(x) - 2) * 10000

# Que está bem próximo do valor encontrado.

#----Questão 7----

# Comparação do valor absoluto das estatísticas de teste obtidas com
# os quantis (1-0.05/2)% da distribuição t com 8 graus de liberdade:

sum((abs(estat_teste_beta_1_2) > qt(1 - 0.05/2, length(x) - 2)) |
    (abs(estat_teste_beta_0) > qt(1 - 0.05/2, length(x) - 2)))

# Aqui temos 716 simulações em que pelo menos uma das duas hipóteses foi rejeitada.
