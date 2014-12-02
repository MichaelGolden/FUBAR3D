# 0 ≤ β ≤ 5
# phase transition at β ≈ 2.5

# run ising model
function isingSim(S::SparseMatrixCSC,β::Float64)
    # run ising model given structure (weighted graph)

    # return vector of states
end

# β repitition & selection
function repβ(S::SparseMatrixCSC,β::Float16=1:.1:5,rep::Int32=10)
#  rep>typemax(typeof(rep)) && error("rep must be less than",typemax(typeof(rep)))
  Mβ = Array{Int8,3}
  Mβ = [i for i=isingSim(S,β),r=1:rep,βi=β]
end
