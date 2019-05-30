library(data.table)
library(ggplot2)

# for testing purposes
# source("/da/dmp/cb/wegmare1/scRNASeq/workflow_devel/scRNASeq_workflow/code/scRNASeq_pipeline_functions.R")
# plot_dt = data.table(tSNE1 = c(1,2,3,4), tSNE2=c(1,2,3,4), col1 = c(1,1,1.5,2), col2 = c(1,1,1.5,2),factor = factor(c(1,1,2,2)))

ui <- fluidPage(

  titlePanel("Marker gene expression"),
  
  sidebarLayout(
    
    sidebarPanel(

      selectInput("color", label = "Choose marker:", 
                  choices = as.list(names(plot_dt)[-c(1,2)]),
                  selected = names(plot_dt)[3]),
    
      sliderInput("size", label = "Point size:",
                min = 1, max = 10, value = 3),
      
      sliderInput("alpha", label = "Transparency:",
                  min = 0.1, max = 1, value = 0.8),
      
      selectInput("pal", label = "Color palette:", 
                  choices = list('default', 'red-blue', 'heat'),
                  selected = 'default')
      
    ),

    mainPanel(
      h3(textOutput("selected_col")),
      plotOutput("plot"))
  )
  
)  
 

# Define server logic ----
server <- function(input, output) {
  
  output$selected_col <- renderText({ 
    paste("tSNE map, colored by", input$color)
  })
  
  
  output$plot  = renderPlot({
    switch(input$pal,
    'default' = generic_scatterplot(plot_dt,
                        x_col = "tSNE1",
                        y_col="tSNE2",
                        color = input$color,
                        abs_size = input$size,
                        alpha = input$alpha),
    
    'red-blue' = generic_scatterplot(plot_dt,
                          x_col = "tSNE1",
                          y_col="tSNE2",
                          color = input$color,
                          abs_size = input$size,
                          alpha = input$alpha) + scale_color_distiller(palette='RdBu'),
    
    'heat' = generic_scatterplot(plot_dt,
                                     x_col = "tSNE1",
                                     y_col="tSNE2",
                                     color = input$color,
                                     abs_size = input$size,
                                     alpha = input$alpha) + scale_color_distiller(palette='RdYlBu')
    )
    
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)