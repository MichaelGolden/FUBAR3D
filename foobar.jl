using Distributions

function likelihood(A::Int8, mat::Array{Float64}, G::MvNormal)   #::Float64
    sum = 0
    dim = int(sqrt(length(mat)))
    for i in length(mat)
        dn = i / dim
        ds = i % dim
        sum += mat[i]*pdf(G,[dN, dS])
    end
    return sum
end

function logLikelihood(A::Array{Int8}, mat::Array{Float64}, G::MvNormal)   #::Float64
    sum = 0
    for Ai in A
        sum += log(likelihood(Ai, mat, G))
    end
end

function discretize(G::MvNormal, inputSpace::Array{Float64})   #::MvNormal
    dict = {}
    normalization = 0
    for input in inputSpace
        pdf = pdf(G,input)
        dict[input] = pdf
        normalization += pdf
    end
    normalization = 1.0 / normalization

    total = 0
    for input, pdf in dict
        pdf = pdf / normalization
        dict[input] = pdf
        total += pdf
    end

    assert total == 1

    # convert dict to a distribution object, perhaps
    return dict
end



function update()
  function likelihood(A::Int8, mi::Float64, d::MvNormal)   #::Float64
 #   sum = 0
 #   dim = int(sqrt(length(mat)))
 #   for i in length(mat)
 #       dn = i / dim
 #       ds = i % dim
#        sum += mat[i]*pdf(G,[dN, dS])
        sum = Ai*pdf(d,mi)
 #   end
  #  return sum
end

# A: Array{Int8,1}, 400x1
# mat: Array{Float64,2}, 400,2
#
function logLikelihood(A::Array{Float64,1}, mat::Array{Float64,2}, d::MvNormal)   #::Float64
    sum = 0
    for i in 1:length(A)
        sum += log(likelihood(A[i], mat[i,], d))
    end
end
  function test_lik()
  d = MultiVnorm([2.5,2.5],[1.,1.],.1)
  gridInfoFile = "datasets/rhodopsin.nwk.grid_info"
  h = getHeights(d,gridInfoFile)
  include("FUBARDataset.jl")
  grid = loadDataset(gridInfoFile)[3]
  println(typeof(h))
  L = logLikelihood(h,grid,d)
  println(L)
end

test_lik()
end
