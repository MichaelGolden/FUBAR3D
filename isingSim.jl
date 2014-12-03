include("MCMC.jl")
# 0 ≤ β ≤ 5
# phase transition at β ≈ 2.5

# run ising model
function isingSim(S::SparseMatrixCSC{Float64,Int64},β::Float64,maxIter=1e4)
    # run ising model given structure (weighted graph)
    mainMCMC(S,β,maxIter)
    # return vector of states  Array{Int8,1}
end

# β repitition & selection
function repβ(S::SparseMatrixCSC{Float64,Int64},β::FloatRange{Float64}=1:.1:5,rep::Int64=10)
#  rep>typemax(typeof(rep)) && error("rep must be less than",typemax(typeof(rep)))
  vβ = Array{Int8,1}
  for b in β
    for r in 1:rep
#      for i in isingSim(S,b)
        push!(vβ,isingSim(S,b))
 #     end
    end
  end
  Mβ = reshape(vβ,(size(S)[1],rep,length(β)))
#  Mβ = sizehint(Array(Int8,1), size(S)[1]*rep*length(β))
#  Mβ = [i for i=isingSim(S,b),r=1:rep,b=β]
end

### test
out = repβ(sparse(readdlm("datasets/distance.matrix.thresh15.csv",' ')))

