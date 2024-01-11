library(shiny)
library(smam)
library(tidyverse)
library(shinycssloaders)
library(gganimate)
library(gifski)
library(av)

# devtools::install_version("MASS", "7.3-51.1")



# Define UI for application that draws a histogram
ui <- navbarPage(title = "smam web-app",
    # mrme
    tabPanel(title = "Moving-resting process with measurement error",
            fluidPage(

              fluidRow(
                column(6, code("Fit a Moving-Resting Model with Measurement Error. The measurement error is modeled by Guassian noise."))
              ),
              tags$hr(),
              sidebarLayout(
                sidebarPanel(
                  fileInput("file_mrme", HTML("select animal location file<br/>(A .csv file including ONLY 3 NAMED columns with first column is <tt>time</tt> (hr) and rest two columns are <tt>north</tt>, <tt>east</tt> coordinator (km))"),
                            multiple = FALSE,
                            accept = ".csv"),
                  numericInput("start_mrme_lamM", "starting value of lamM", 5,
                               min = 0),
                  numericInput("start_mrme_lamR", "starting value of lamR", 1,
                               min = 0),
                  numericInput("start_mrme_sigma", "starting value of sigma", 1,
                               min = 0),
                  numericInput("start_mrme_sig_err", "starting value of sig_err", 0.01,
                               min = 0),
                  actionButton("run_mrme", "fit MRME")
                ),
                mainPanel(
                  h2("GPS coordinates of animal"),
                  plotOutput('plot_mrme'),
                  
                  tags$hr(),
                  
                  h2("Trace plot of animal"),
                  imageOutput("plot_mrme2"),
                  
                  tags$hr(),
                  
                  
                  h2("Display input dataset"),
                  DT::dataTableOutput("table_mrme"),
                  
                  tags$hr(),
                  
                  h2("MRME process estimation result"),
                  textOutput('text_mrme') %>% withSpinner(color="#0dc5c1")
                  )
              )
            )),
    # mm
    tabPanel(title = "Moving-moving process",
             fluidPage(
               fluidRow(
                 column(6, code("Fit a Moving-Moving Model with 2
                                Embedded Brownian Motion with animal
                                movement data at discretely observation
                                times by maximizing a full likelihood
                                constructed from the marginal density
                                of increment."))
               ),
               tags$hr(),
               sidebarLayout(
                 sidebarPanel(
                   fileInput("file_mm", HTML("select animal location file<br/>(A .csv file including ONLY 3 NAMED columns with first column is <tt>time</tt> (hr) and rest two columns are <tt>north</tt>, <tt>east</tt> coordinator (km))"),
                             multiple = FALSE,
                             accept = ".csv"),
                   numericInput("start_mm_lamM1", "starting value of lamM1", 1,
                                min = 0),
                   numericInput("start_mm_lamM2", "starting value of lamM2", 1,
                                min = 0),
                   numericInput("start_mm_sigma1", "starting value of sigma1", 1,
                                min = 0),
                   numericInput("start_mm_sigma2", "starting value of sigma2", 0.1,
                                min = 0),
                   actionButton("run_mm", "fit MM")
                 ),
                 mainPanel(
                   h2("GPS coordinates of animal"),
                   plotOutput('plot_mm'),
                   
                   tags$hr(),
                   
                   h2("Display input dataset"),
                   DT::dataTableOutput("table_mm"),
                   
                   tags$hr(),
                   
                   h2("MM process estimation result"),
                   textOutput('text_mm') %>% withSpinner(color="#0dc5c1")
                 )
               )
             ))
)



server <- function(input, output) {

  ## MRME
  input_data_mrme <- reactive({
    if (is.null(input$file_mrme)) return(NULL)
    read.csv(file = input$file_mrme$datapath)
  })
  
  output$table_mrme <- DT::renderDataTable({
    
    # render only if there is data available
    req(input_data_mrme())
    
    # reactives are only callable inside an reactive context like render
    data <- input_data_mrme()
    print(data)
    data
  })
  
  output$plot_mrme <- renderPlot({
    req(input_data_mrme())
    data <- input_data_mrme()
    
    data <- rbind(c(0, 0, 0),
                  apply(data, 2, function(x) cumsum(diff(x))))
    data <- as.data.frame(data)
    data %>%
      ggplot() +
      geom_line(aes_string(x = colnames(data)[1], y = colnames(data)[2])) +
      geom_line(aes_string(x = colnames(data)[1], y = colnames(data)[3])) +
      theme_bw() +
      xlab("time") +
      ylab("coordinates")
  })
  
  output$plot_mrme2 <- renderImage({
    req(input_data_mrme())
    data <- input_data_mrme()
    
    data <- rbind(c(0, 0, 0),
                  apply(data, 2, function(x) cumsum(diff(x))))
    data <- as.data.frame(data)
    
    ani_plot_mrme2 <- ggplot(data, aes(x = east, y = north)) +
      geom_point(size = 3) +
      theme_bw() + coord_fixed(ratio = 1) +
      labs(title = 'time/hr.: {frame_time}', x = 'east-west/km', y = 'north-south/km') +
      transition_time(time) +
      shadow_mark(colour = 'grey') +
      ease_aes('linear')
    
    anim_save("ani_plot_mrme2.gif", animate(ani_plot_mrme2, renderer = gifski_renderer()))
    
    list(src = "ani_plot_mrme2.gif",
         contentType = 'image/gif'
         # width = 400,
         # height = 300,
         # alt = "This is alternate text"
         )
    
  }, deleteFile = TRUE)
  
  dofit_mrme <- eventReactive(input$run_mrme, {
    req(input_data_mrme())
    data <- input_data_mrme()

    result <- fitMRME(data,
                      start = c(input$start_mrme_lamM,
                                input$start_mrme_lamR,
                                input$start_mrme_sigma,
                                input$start_mrme_sig_err))

    result
  })
  


  output$text_mrme <- renderText({
    result <- dofit_mrme()

    paste0("Estimation result is ", paste(format(result$estimate, digits = 5), collapse = " "),
           ". Log likelihood is ", format(result$loglik, digits = 5),
           ". Convergence status is ", result$convergence, ".")
  })

  
  
  ## MM
  input_data_mm <- reactive({
    if (is.null(input$file_mm)) return(NULL)
    read.csv(file = input$file_mm$datapath)
  })
  
  output$table_mm <- DT::renderDataTable({
    
    # render only if there is data available
    req(input_data_mm())
    
    # reactives are only callable inside an reactive context like render
    data <- input_data_mm()
    print(data)
    data
  })
  
  output$plot_mm <- renderPlot({
    req(input_data_mm())
    data <- input_data_mm()
    
    data <- rbind(c(0, 0, 0),
                  apply(data, 2, function(x) cumsum(diff(x))))
    data <- as.data.frame(data)
    data %>%
      ggplot() +
      geom_line(aes_string(x = colnames(data)[1], y = colnames(data)[2])) +
      geom_line(aes_string(x = colnames(data)[1], y = colnames(data)[3])) +
      theme_bw() +
      xlab("time") +
      ylab("coordinates")
  })
  
  
  dofit_mm <- eventReactive(input$run_mm, {
    req(input_data_mm())
    data <- input_data_mm()
    
    result <- fitMM(data,
                    start = c(input$start_mm_lamM1,
                              input$start_mm_lamM2,
                              input$start_mm_sigma1,
                              input$start_mm_sigma2))
    
    result
  })
  
  
  output$text_mm <- renderText({
    result <- dofit_mm()
    
    paste0("Estimation result is ", paste(format(result$estimate, digits = 5), collapse = " "),
           ". Log likelihood is ", format(result$loglik, digits = 5),
           ". Convergence status is ", result$convergence, ".")
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
