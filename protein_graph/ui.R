library(shiny)
library(shinyRGL)

# Define UI for application that draws a histogram
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Shiny WebGL!"),
  
  # Sidebar with a slider input for number of points
  sidebarPanel(
    sliderInput("pts", 
                "Number of points:", 
                min = 10, 
                max = 1000, 
                value = 250),
    HTML("<hr />"),
    helpText(HTML("Created using <a href = \"http://github.com/trestletech/shinyRGL\">shinyRGL</a>."))
  ),
  
  # Show the generated 3d scatterplot
  mainPanel(
    webGLOutput("sctPlot")
  )
))