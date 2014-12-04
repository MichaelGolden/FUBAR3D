include("MCMC.jl")
# 0 ≤ β ≤ 5
# phase transition at β ≈ 2.5

# run ising model
function isingSim(S::SparseMatrixCSC{Float64,Int64},β::Float64,grid,maxIter=1e4,repNum=0)
    # run ising model given structure (weighted graph)
    mainMCMC(S,grid,β,maxIter,repNum)
    # return vector of states  Array{Int8,1}
end

# β repitition & selection
function repβ(S::SparseMatrixCSC{Float64,Int64},grid,β::FloatRange{Float64}=1:.1:5,rep::Int64=10)
  Mβ = zeros(Int,(size(S)[1],rep,length(β)))
  for b in 1:length(β)
    for r in 1:rep
      i = isingSim(S,β[b],grid,1e4,r)
      Mβ[:,r,b] += i
    end
  end
  Mβ
end

### test
function test()
  thresh = 5.1
  R = readdlm("datasets/distance.matrix.thresh65.csv",' ')
  for i in 1:size(R)[1]
    for j in 1:size(R)[2]
      if R[i,j] > thresh
        R[i,j] = 0
      end
    end
  end
  S = sparse(R)
  β = 0:.5:3
  out = repβ(S,"datasets/rhodopsin.nwk.grid_info",β,5)
  writedlm3D(out,"out/test.out",β,',')
end

function writedlm3D(m,title,d3Titles,dlm)
  if length(d3Titles) != size(m)[3]
    error("length(d3Titles) != size(m)[3] :: d3Titles is the titles of the 3rd dimension of the array.")
  end
  for i in 1:(size(m)[3])
    writedlm(string(title,".",d3Titles[i],".csv"),m[:,:,i],dlm)
  end
end

