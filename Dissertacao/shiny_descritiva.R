library(shiny)
library(leaflet)
library(ggplot2)
library(dplyr)
library(sf)
library(terra)
library(spatstat)
library(lubridate)
library(parallel)

#----Funções----

ui <- fluidPage(
    titlePanel("Dashboard ETAS: Simulação, Estimação e Descendência"),
    sidebarLayout(
        sidebarPanel(
            tabsetPanel(
                tabPanel("Simulação",
                         numericInput("seed", "Seed para simulação", 5050),
                         h4("Parâmetros do Modelo (Theta)"),
                         numericInput("mu", "mu (Background)", 0.01, step = 0.01),
                         numericInput("A", "A (Produtividade)", 0.5, step = 0.1),
                         numericInput("alpha", "alpha (Eficiência)", 1.0, step = 0.1),
                         numericInput("c_par", "c (Tempo)", 0.2, step = 0.05),
                         numericInput("p_par", "p (Omori)", 1.15, step = 0.05),
                         numericInput("d_par", "D (Espaço)", 0.02, step = 0.01),
                         sliderInput("t_max", "Tempo Máximo (Anos)", 2, 50, 10),
                         actionButton("simulate_btn", "Simular Base de Dados", class = "btn-primary")
                         ),
                tabPanel("Filtros",
                         sliderInput("time_range", "Janela Temporal:", 0, 10, c(0, 10)),
                         sliderInput("mag_filter", "Magnitude Mínima:", 4, 9, 4, step = 0.5),
                         selectizeInput("search_id", "Buscar Terremoto por ID:",
                                        choices = NULL,
                                        options = list(placeholder = 'Digite o ID do evento'))
                         )
            ),
            hr(),
            actionButton("estimate_btn", "Estimar Background (Zhuang)", class = "btn-warning"),
            helpText("O clique no mapa filtra a descendência do terremoto selecionado.")
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("Mapa & Descendência",
                         leafletOutput("map_plot", height = 500),
                         wellPanel(textOutput("click_info"))),
                tabPanel("Tabela",
                         tableOutput("tab_rank")),
                tabPanel("Série Temporal",
                         plotOutput("timeline_plot")),
                tabPanel("Estimação",
                         plotOutput("raster_plot"),
                         plotOutput("prob_plot"))
            )
        )
    )
)

server <- function(input, output, session) {
    # Valores para armazenar os dados:
    v <- reactiveValues(
        raw_data = NULL,
        filtered_data = NULL,
        estimation_raster = NULL,
        selected_id = NULL,
        descendants = NULL
    )
    # Simulação dos eventos:
    observeEvent(input$simulate_btn, {
        withProgress(message = 'Simulando ETAS...', {
            # Seed definida pelo usuário:
            set.seed(input$seed)
            sim <- simulate_hawkes_etas_zhuang(
                T_max = input$t_max,
                mu = input$mu,
                A = input$A,
                u_fn = rast_peru,
                u_max = max(im_peru),
                alpha = input$alpha,
                c_par = input$c_par,
                p_par = input$p_par,
                d_par = input$d_par,
                beta = 2, m_min = 4,
                # Limites baseados no Peru ou genéricos
                X_lim = c(min(peru.quakes$longlat.coord$long), max(peru.quakes$longlat.coord$long)),
                Y_lim = c(min(peru.quakes$longlat.coord$lat), max(peru.quakes$longlat.coord$lat))
                )
            if(!is.null(sim)) {
                sim$id <- rownames(sim)
                v$raw_data <- sim
                updateSelectizeInput(session, "search_id",
                                     choices = c("", sim$id),
                                     selected = "",
                                     server = TRUE)
                updateSliderInput(session, "time_range", max = input$t_max, value = c(0, input$t_max))
            }
        })
    })
    # Procura de eventos por ID:
    observeEvent(input$search_id, {
        if (input$search_id == "") {
            v$selected_id <- NULL
            v$descendants <- NULL
        } else {
            v$selected_id <- as.character(input$search_id)
            find_descendants <- function(data, parent_id) {
                children <- data %>% filter(as.character(parent) == as.character(parent_id))
                return(children$id)
            }
            v$descendants <- find_descendants(v$raw_data, v$selected_id)
        }
    })
    # Filtro temporal:
    observe({
        req(v$raw_data)
        v$filtered_data <- v$raw_data %>%
            filter(t >= input$time_range[1], t <= input$time_range[2],
                   m >= input$mag_filter)
    })
    # Filtro por clique:
    observeEvent(input$map_plot_marker_click, {
        click <- input$map_plot_marker_click
        req(click$id)
        v$selected_id <- as.character(click$id)
        find_descendants <- function(data, parent_id) {
            children <- data %>% filter(as.character(parent) == as.character(parent_id))
            return(children$id)
        }
        v$descendants <- find_descendants(v$raw_data, v$selected_id)
        updateSelectizeInput(session, "search_id", selected = v$selected_id)
    })
    # Mapa com eventos completos:
    output$map_plot <- renderLeaflet({
        req(v$filtered_data)
        df <- v$filtered_data
        leaflet(df) %>%
            addTiles() %>%
            addCircleMarkers(~x, ~y, layerId = ~id,
                             radius = ~m, color = "blue",
                             stroke = TRUE, fillOpacity = 0.6,
                             popup = ~paste("ID:", id,
                                            "Time:", round(t, 2),
                                            "Mag:", round(m, 2),
                                            "Gen:", gen))
    })
    # Atualização do mapa conforme filtros:
    observe({
        req(v$filtered_data)
        df <- v$filtered_data
        df$color <- "blue"
        df$opacity_val <- 0.6
        if(!is.null(v$selected_id)) {
            df$color <- "grey"
            df$opacity_val <- 0.15
            # Evento Principal (Parent):
            df$color[as.character(df$id) == as.character(v$selected_id)] <- "red"
            df$opacity_val[as.character(df$id) == as.character(v$selected_id)] <- 1.0
            # Descendentes (Offspring):
            df$color[as.character(df$id) %in% as.character(v$descendants)] <- "orange"
            df$opacity_val[as.character(df$id) %in% as.character(v$descendants)] <- 0.9
        }
        # Atualiza o mapa conforme cliques:
        leafletProxy("map_plot", data = df) %>%
            clearMarkers() %>%
            addCircleMarkers(~x, ~y, layerId = ~id,
                             radius = ~m,
                             color = ~color,
                             stroke = TRUE,
                             weight = 1,
                             fillOpacity = ~opacity_val,
                             opacity = ~opacity_val,
                             popup = ~paste("ID:", id,
                                            "Time:", round(t, 2),
                                            "Mag:", round(m, 2),
                                            "Gen:", gen,
                                            "Parent:", parent))
    })
    # Informações de parents e offsprings para o evento selecionado:
    output$click_info <- renderText({
        if(is.null(v$selected_id)) "Clique em um terremoto para ver seus descendentes diretos."
        else paste("Terremoto selecionado:", v$selected_id, "| Descendentes encontrados:", length(v$descendants))
    })
    # Timeline (futuramente trocar para série temporal da intensidade condicional integrada no espaço):
    output$timeline_plot <- renderPlot({
        req(v$filtered_data)
        ggplot(v$filtered_data, aes(x = t, y = m)) +
            geom_segment(aes(xend = t, yend = 4), alpha = 0.3) +
            geom_point(aes(color = as.factor(gen))) +
            labs(title = "Série Temporal (Colorido por Geração)", x = "Tempo", y = "Magnitude", color = "Geração") +
            theme_minimal()
    })
    output$tab_rank <- renderTable({
        req(v$raw_data)
        df <- v$raw_data
        df %>%
            count(parent) %>%
            arrange(desc(n)) %>%
            head(10) %>%
            select("Parents" = parent,
                   "N_Offspring" = n)
    })
    # Ajuste do background rate:
    observeEvent(input$estimate_btn, {
        req(v$raw_data)
        withProgress(message = 'Rodando Algoritmo de Zhuang...', value = 0, {
            df_est <- v$raw_data %>%
                mutate(
                    time = t,
                    mag = m,
                    long = x,
                    lat = y
                )
            catalogo <- catalog(
                df_est %>%
                mutate(t_years = dyears(t),
                       data = POSIXct(1) + t_years,
                       date = date(data),
                       time_fmt = format(data, format = "%H:%M:%S"))
            )
            n_events <- nrow(df_est)
            theta <- c("A" = input$A, "alpha" = input$alpha, "c" = input$c_par,
                       "p" = input$p_par, "D" = input$d_par, "beta" = 2)
            r <- terra::rast(nrows = 128, ncols = 128,
                             xmin = min(df_est$x),
                             xmax = max(df_est$x),
                             ymin = min(df_est$y),
                             ymax = max(df_est$y))
            values(r) <- 1
            d_j_simulado <- d_j(catalogo$longlat.coord[, c("long", "lat")])
            rasters <- list()
            rasters[[1]] <- r
            num_cores <- if(.Platform$OS.type == "windows") 1 else max(1, detectCores() - 2)
            calcular_probabilidades <- function(r_atual) {
                probs <- mclapply(1:n_events, function(j) {
                    if (j == 1) return(0)
                    i_indices <- 1:(j - 1)
                    denom <- lambda_total(
                        df_est$time[j], df_est$x[j], df_est$y[j], df_est$mag[j],
                        4, df_est, theta, r_atual
                    )
                    numeradores <- kappa_func(df_est$mag[i_indices], 4, theta["A"], theta["alpha"]) *
                        g_func(df_est$time[j] - df_est$time[i_indices], theta["c"], theta["p"]) *
                        f_spatial(
                        (df_est$x[j] - df_est$x[i_indices])^2 + (df_est$y[j] - df_est$y[i_indices])^2,
                        df_est$mag[i_indices], 4, theta["D"], theta["alpha"]
                        )
                    return(sum(numeradores / denom))
                }, mc.cores = num_cores)
                return(unlist(probs))
            }
            # -----------------------------------------------------------------
            incProgress(0.2, detail = "Calculando Iteração 1...")
            # Primeira iteração:
            prob_total <- calcular_probabilidades(r)
            rho_background <- pmax(0, 1 - prob_total)
            values(r) <- estimate_background_kernel(
                seq(min(df_est$long), max(df_est$long), length.out = 128),
                seq(min(df_est$lat), max(df_est$lat), length.out = 128),
                data.frame(x = df_est$x, y = df_est$y, rho = rho_background),
                d_j_simulado,
                T_total = max(catalogo$rtperiod)
            )
            r <- flip(r, direction = "vertical")
            rasters[[2]] <- r
            # Loop de Convergência
            i <- 2
            max_iter <- 10
            while (any(abs(values(rasters[[i]]) - values(rasters[[i - 1]])) > 10e-6) && i < max_iter) {
                incProgress(1/max_iter, detail = paste("Calculando Iteração", i, "..."))
                prob_total <- calcular_probabilidades(rasters[[i]])
                rho_background <- pmax(0, 1 - prob_total)
                r_loop <- rasters[[i]]
                values(r_loop) <- estimate_background_kernel(
                    seq(min(df_est$long), max(df_est$long), length.out = 128),
                    seq(min(df_est$lat), max(df_est$lat), length.out = 128),
                    data.frame(x = df_est$x, y = df_est$y, rho = rho_background),
                    d_j_simulado,
                    T_total = max(catalogo$rtperiod)
                )
                r_loop <- flip(r_loop, direction = "vertical")
                rasters[[i + 1]] <- r_loop
                i <- i + 1
            }
            incProgress(1, detail = "Finalizado!")
            v$estimation_raster <- rasters[[i]]
            v$probabilidades_totais <- prob_total
        })
    })
    output$raster_plot <- renderPlot({
        req(v$estimation_raster)
        plot(v$estimation_raster, main = "Intensidade Estimada do Background")
    })
    output$prob_plot <- renderPlot({
        req(v$probabilidades_totais)
        plot(v$probabilidades_totais, main = "Probabilidades Totais por Evento")
    })
}

shinyApp(ui, server)
