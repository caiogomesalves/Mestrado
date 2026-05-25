#----Pacotes---------------------------------
library(tidyverse)

#----Simulação-------------------------------

# Seed para reprodutibilidade:
set.seed(5050)

n <- 100
datas <- as.Date(c("1960-01-01", "2015-12-31"))

pacientes <- tibble(
    id_paciente = paste0("P",
                         sprintf("%03d", 1:n)),
    sexo = factor(sample(c("M", "F", NA), n, T,
                         c(0.45, 0.45, 0.1)),
                  levels = c("M", "F")),
    data_nascimento = sample(as.Date(datas[1]:datas[2]), n,
                             replace = F),
    peso_kg = case_when(
        sexo == "M" ~ rnorm(n, 80, 5),
        sexo == "F" ~ rnorm(n, 55, 4),
        .default = rnorm(n, 70, 2)
    )
)

n_exames <- 300
datas2 <- as.Date(c("2024-01-01", "2024-12-31"))

exames <- tibble(
    id_paciente = sample(pacientes$id_paciente, n_exames, T),
    data_exame = sample(as.Date(datas2[1]:datas2[2]), n_exames, F),
    unidade = sample(c("mg/Dl", "mmol/L"), n_exames, T),

    )
