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
peptides <- lapply(s, function(x) points[,c(x, x + 1, x + 2)])

# 3-d plot the coordinates at each site for each peptide
plot_peptide <- function(p)
{
  plot3d(p[,1], p[,2], p[,3], col=sample(colors(), 1), 
         add=TRUE, size=5)
  ts.surf <- t(convhulln(p))
  rgl.triangles(p[ts.surf,1],p[ts.surf,2],p[ts.surf,3],
               col="blue",alpha=0.05)
}
lapply(peptides, plot_peptide)

max_dist = 10

# function to draw segments between all pairs of close points
connect_peptides <- function(peptide1, peptide2)
{
  #g1 = g[[1]]
  #g2 = g[[2]]
  
  probs = replicate(nrow(peptide1), runif(nrow(peptide1)))
  close = fields.rdist.near(peptide1, peptide2, delta=max_dist, 
                            max.points=nrow(peptide1)^2)
  close_inds = close$ind
  close_inds <- close_inds[close_inds[,1] != close_inds[,2],]
  
  close1 = close_inds[,1]
  close2 = close_inds[,2]
  
  peptide1 = peptide1[close1,]
  peptide2 = peptide2[close2,]
  
  
  probs = probs[close1, close2]
  probs = diag(probs)
  ncolors = 20
  ramp <- heat.colors(ncolors)
  bins <- seq(0, 1, 1/ncolors)
  colors <- .bincode(probs, bins, include.lowest=TRUE)
  colors <- ramp[colors]
  
  
  segments3d(c(peptide1[,1], peptide2[,1]), 
             c(peptide1[,2], peptide2[,2]), 
             c(peptide1[,3], peptide2[,3]), 
             add=TRUE, col=colors)
}

# iterate through every combination of two peptides
# for each combination apply the function to connect close sites
# between those two peptide pairs
for (peptide1 in peptides)
{
  for (peptide2 in peptides)
  {
    connect_peptides(peptide1, peptide2)
  }
}
