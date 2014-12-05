library(rgl)
library(fields)
library(pdist)
library(magic)
library(abind)
library(geometry)
#library(combinat)
library(spam)
library(grid)
library(grDevices)

points <- read.table("datasets/hiv1.coords", sep="\t", header=TRUE)
points <- na.omit(points)

s <- seq(1, ncol(points), by=3)
# split the 3n columns into n data frames
# data frame columns contain x, y, z coordinates for the peptide
# rows contain the xyz coordinate for each site
frames <- lapply(s, function(x) points[,c(x, x + 1, x + 2)])

# 3-d plot the coordinates at each site for each peptide
plot_peptide <- function(p)
{
  plot3d(p[,1], p[,2], p[,3], col=sample(colors(), 1), add=TRUE)
  ts.surf <- t(convhulln(p))
  #rgl.triangles(p[ts.surf,1],p[ts.surf,2],p[ts.surf,3],
   #            col="blue",alpha=0.01)
}
lapply(frames, plot_peptide)

max_dist = 10

# function to draw segments between all pairs of close points
connect_groups <- function(g)
{
  g1 = g[[1]]
  g2 = g[[2]]
  
  colors = replicate(nrow(g1), runif(nrow(g1)))
  close = fields.rdist.near(g1, g2, delta=max_dist, 
                            max.points=nrow(g1)^2)
  close_inds = close$ind
  close_inds <- close_inds[close_inds[,1] != close_inds[,2],]
  
  close1 = close_inds[,1]
  close2 = close_inds[,2]
  
  g1 = g1[close1,]
  g2 = g2[close2,]
  
  colors = colors[close1, close2]
  colors = diag(colors)
  ncolors = 10000
  ramp <- heat.colors(ncolors)
  bins <- seq(0, 1, 1/ncolors)
  colors <- .bincode(colors, bins, include.lowest=TRUE)
  colors <- ramp[colors]
  
  
  segments3d(c(g1[,1], g2[,1]), c(g1[,2], g2[,2]), c(g1[,3], g2[,3]), 
             add=TRUE, col=colors)
}

# iterate through every combination of two peptides
# for each combination apply the function to connect close sites
# between those two peptide pairs
for (frame1 in frames)
{
  for (frame2 in frames)
  {
    connect_groups(list(frame1, frame2))
  }
}
