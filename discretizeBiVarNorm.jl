using Plotly

#metropolis hastings moves


# Discretize given distribution onto DnDs space
## input: Distribution, grid
## output: Normalized heights
function discMVNorm(d::MultiVnorm,grid::Array{Float64,2},rankgrid::Bool=true)
  if rankgrid
    grid = RankGrid(grid)
  end
  # iterate through grid
    # evaluate height (probability | Distribution)
  h = [pdf(d,[grid[i,1],grid[i,2]]) for i=1:size(grid)[1]]
  # Normalize so all the heights sum to 1
  Σh = sum(h)
  hNormalized = [float64(h[i]/Σh) for i=1:length(h)]
end

#rank
function RankGrid(g::Array{Float64,2})
  out = zeros(size(g))
  r = 1:length(unique(g[:,1]))
  l = length(r)
  for i in r
    out[(l*(i-1)+1):(l*i),1] = i
    out[(l*(i-1)+1):(l*i),2] = r
  end
  out
end

function test_rank()
    d = MultiVnorm([2.5,2.5],[1.,1.],.1)
  include("FUBARDataset.jl")
  grid = loadDataset("datasets/rhodopsin.nwk.grid_info")[3]
    h = discMVNorm(d,grid)
end

a=test_rank()
typeof(a)
println(a')
#cd("/Users/benjaminkellman/Documents/Grad/Labs/Pond/WeekPaper/FUBAR3D/")
#println(test_rank())
