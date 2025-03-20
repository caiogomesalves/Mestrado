# Questão 5

## Número de mulheres
n_f <- 241945

## Número de homens
n_m <- 251527

## Teste exato de fischer (para fins de comparação):
binom.test(c(n_f, n_m), p = 0.5, alternative = "greater")

## Considerando que os dados sejam binomiais, a distribuição tem a seguinte forma:
plot(seq(0, 1, length.out = 1000), dbeta(seq(0, 1, length.out = 1000), 241916, 251527), type = "l")

## Para fazermos inferência bayesiana, consideremos uma priori U(0, 1) que é exatamente uma Beta(1, 1):
plot(dbeta(seq(0, 1, length.out = 1000), 1, 1), type = "l")

## Sabemos que a posteriori será conjugada, sendo Beta(n_f + 1, n_m + 1).
## Assim, para fazermos inferência sobre essa distribuição por meio de simulação, devemos inicialmente amostrar dessa Beta.

## Seed para gerar amostras da posteriori:
set.seed(5050)

## Gera um milhão de amostras da Beta
amostra_theta <- rbeta(1e6, n_f + 1, n_m + 1)

## Proporção de valores que são maiores do que 0.5:
mean(amostra_theta > 0.5)

## O formato da posteriori é a seguinte:
plot(density(amostra_theta))

## Note o eixo x no gráfico anterior. Uma representação mais correta seria:
plot(density(amostra_theta), xlim = c(0, 1))
