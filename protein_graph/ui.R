library(shiny)
library(shinyRGL)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Show the generated 3d scatterplot
  mainPanel(
    webGLOutput("sctPlot")
  )
))