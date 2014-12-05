library(rgl)
library(fields)
library(pdist)

points <- read.table("datasets/hiv1.coords", sep="\t", header=TRUE)
points <- na.omit(points)

s <- seq(1, ncol(points), by=3)
# split the 3n columns into n data frames
# data frame columns contain 
frames <- lapply(s, function(x) points[,c(x, x + 1, x + 2)])

lapply(frames, function(f) 
  {plot3d(f[,1], f[,2], f[,3], col=sample(colors(), 1), add=TRUE)})

max_dist = 20

connect_groups <- function(g)
{
  g1 = g[[1]]
  g2 = g[[2]]
  
  #print(g1)
  #print(g2)
  
  close = fields.rdist.near(g1, g2, delta=max_dist, 
                            max.points=nrow(g1)^2)
  close_inds = close$ind
  #print(close_inds)
  close_inds <- close_inds[close_inds[,1] != close_inds[,2],]
  
  close1 = close_inds[,1]
  close2 = close_inds[,2]
  
  #print(close1)
  #print(close2)
  
  g1 = g1[close1,]
  g2 = g2[close2,]
  
  #print(g1)
  #print(g1)
  
  #print(g1[,3])
  #print(g2[,3])
  
  colors = replicate(nrow(g1), runif(nrow(g1)))
  colors = colors[close1, close2]
  colors = diag(colors)
  
  #print(length(colors))
  #print(nrow(g1))
  
  segments3d(c(g1[,1], g2[,1]), c(g1[,2], g2[,2]), c(g1[,3], g2[,3]), 
             add=TRUE, col=10*colors)
  return(1)
}

combs <- combn(frames, 2, FUN=connect_groups)

