
pacman::p_load("renv", "here", "tidyverse","shiny","brms","patchwork")

# Load your utility functions
source(here::here("Shiny","utility_shiny.R"))

# get contngency space from one participant (two options of contingency1.csv or contingency2.csv)
u = read.csv(here::here("Shiny","Data","contingency1.csv"))
u = u$x

trial = length(u)

# Define UI for the application
ui <- fluidPage(
  
  # Application title
  titlePanel("Shiny App for the various learning models"),
  
  # Define a sidebar layout with input and output
  sidebarLayout(
    sidebarPanel(
      # Navigation bar for different pages
      navlistPanel(
        id = "nav",  # Assign an ID to track the selected tab
        "Learning models",
        tabPanel("hgf", 
                 sliderInput("e_0_hgf", "e_0", min = 0, max = 1, value = 0.5),
                 sliderInput("zeta_hgf", "zeta", min = 0, max = 100, value = 10),
                 sliderInput("omega_hgf", "omega", min = 0, max = 1.5, value = 0.5),
                 sliderInput("pi2_0_hgf", "pi2_0", min = 0, max = 100, value = 0.5)
        ),
        tabPanel("rw", 
                 sliderInput("e_0_rw", "e_0", min = 0, max = 1, value = 0.5),
                 sliderInput("zeta_rw", "zeta", min = 0, max = 100, value = 10),
                 sliderInput("alpha_rw", "alpha", min = 0, max = 1, value = 0.1)
        ),
        tabPanel("ph", 
                 sliderInput("e_0_ph", "e_0", min = 0, max = 1, value = 0.5),
                 sliderInput("zeta_ph", "zeta", min = 0, max = 100, value = 10),
                 sliderInput("S_ph", "S", min = 0, max = 1, value = 0.1),
                 sliderInput("a_0_ph", "a_0", min = 0, max = 1, value = 0.01)
        ),
        tabPanel("su1", 
                 sliderInput("e_0_su1", "e_0", min = 0, max = 1, value = 0.5),
                 sliderInput("zeta_su1", "zeta", min = 0, max = 100, value = 10),
                 sliderInput("mu_su1", "mu", min = 0, max = 15, value = 0.1),
                 sliderInput("Rhat_su1", "Rhat", min = 0, max = 1, value = 0.01),
                 sliderInput("h_0_su1", "h_0", min = -10, max = 10, value = 0.1)
        )
      )
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      uiOutput("main_plot")
    )
  )
)


server <- function(input, output, session) {
  
  output$main_plot <- renderUI({
    tab_name <- input$nav
    plotOutput(tab_name)
  })
  
  output$hgf <- renderPlot({
    data <- get_hgf(u, input$e_0_hgf, input$zeta_hgf, input$omega_hgf, input$pi2_0_hgf)
    
    make_plot(data,"HGF")
  })
  
  output$rw <- renderPlot({
    data <- get_rw(u, input$e_0_rw, input$zeta_rw, input$alpha_rw)
    
    make_plot(data,"RW")
    
  })
  
  output$ph <- renderPlot({
    data <- get_ph(u, input$e_0_ph, input$zeta_ph, input$S_ph,input$a_0_ph)
    
    make_plot(data,"ph")
    
  })
  
  
  output$su1 <- renderPlot({
    data <- get_su1(u, input$e_0_su1, input$zeta_su1, input$mu_su1, input$Rhat_su1, input$h_0_su1)
    
    make_plot(data,"SU1")
  })
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
