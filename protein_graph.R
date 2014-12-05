library(rgl)
library(fields)
library(pdist)

points <- read.table("datasets/hiv1.coords", sep="\t", header=TRUE)
points <- na.omit(points)

s <- seq(1, ncol(points), by=3)
frames <- lapply(s, function(x) points[,c(x, x + 1, x + 2)])

lapply(frames, function(f) 
  {plot3d(f[,1], f[,2], f[,3], col=sample(colors(), 1), add=TRUE)})



delta = 15



