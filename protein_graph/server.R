# Must be executed BEFORE rgl is loaded on headless devices.
options(rgl.useNULL=TRUE)

library(shiny)
library(shinyRGL)
library(rgl)

xPts <- runif(1000, 0, 1)
yPts <- runif(1000, 0, 1)
zPts <- runif(1000, 0, 1)

#' Define server logic required to generate and 3d scatterplot
#' @author Jeff Allen \email{jeff@@trestletech.com}
shinyServer(function(input, output) {
  
  # Expression that generates a rgl scene with a number of points corresponding
  # to the value currently set in the slider.
  output$sctPlot <- renderWebGL({
    #points3d(xPts[1:input$pts],
     #        yPts[1:input$pts],
      #       zPts[1:input$pts])
    #axes3d()
    
    points <- read.table("../datasets/hiv1.coords", sep="\t", header=TRUE)
    points <- na.omit(points)
    npoints <- nrow(points)
    probs <- read.table("../datasets/rhodopsin.nwk.grid_info.sparse.neighbours",
                        sep='\t', header=FALSE)
    probs <- sparseMatrix(i=probs[,1], j=probs[,2], x=probs[,3],
                          dims=c(npoints, npoints))
    
    s <- seq(1, ncol(points), by=3)
    # split the 3n columns into n data frames
    # data frame columns contain x, y, z coordinates for the peptide
    # rows contain the xyz coordinate for each site
    peptides <- lapply(s, function(x) points[,c(x, x + 1, x + 2)])
    
    # 3-d plot the coordinates at each site for each peptide
    plot_peptide <- function(p)
    {
      plot3d(p[,1], p[,2], p[,3], col=sample(colors(), 1), 
             add=TRUE, size=10)
      
      # compute and plot the convex hull for the protein coordinates
      # makes it easier to see links between peptide groups
      # slows down plot-rendering
      ts.surf <- t(convhulln(p))
      #rgl.triangles(p[ts.surf,1],p[ts.surf,2],p[ts.surf,3],
      #            col="blue",alpha=0.05)
    }
    lapply(peptides, plot_peptide)
    
    max_dist = 10
    
    # function to draw segments between all pairs of close points
    connect_peptides <- function(peptide1, peptide2)
    {
      # computes all pairs of points within delta of each other
      close = fields.rdist.near(peptide1, peptide2, delta=max_dist, 
                                max.points=nrow(peptide1)^2)
      close_inds = close$ind
      # remove points from the same row in different peptides
      close_inds <- close_inds[close_inds[,1] != close_inds[,2],]
      
      close1 = close_inds[,1]
      close2 = close_inds[,2]
      
      # get the multiset of points frome each peptide that are
      # close to each other
      peptide1 = peptide1[close1,]
      peptide2 = peptide2[close2,]
      
      # obtain probability of same site between each pair of close
      #peptides
      probs = probs[close1, close2]
      probs = diag(probs)
      
      # create discrete color space, and draw edge colors based on
      # location of the edge probability in color gradient
      #   ncolors = 20
      #   ramp <- heat.colors(ncolors)
      #   bins <- seq(0, 1, 1/ncolors)
      #   colors <- .bincode(probs, bins, include.lowest=TRUE)
      #   colors <- ramp[colors]
      
      # converts the probabilities to a value on a red/black color scale
      # the gb values are set to 0, and r is the probability
      # high probabilities lead to red-ish colors, low probabilities
      # are black
      colors <- rgb(probs, 0, 0)
      #colors <- rgb(1, 1 - probs, 1 -probs)
      
      segments3d(c(peptide1[,1], peptide2[,1]), 
                 c(peptide1[,2], peptide2[,2]), 
                 c(peptide1[,3], peptide2[,3]), 
                 add=TRUE, col=colors, lwd=2)
    }
    
    # connect each peptide pair (including with itself)
    for (peptide1 in peptides)
    {
      for (peptide2 in peptides)
      {
        connect_peptides(peptide1, peptide2)
      }
    }
  })
})