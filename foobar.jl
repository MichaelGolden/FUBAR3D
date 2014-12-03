using Distributions

function likelihood(A::Int8, mat::Array{Float64}, G::Distribution)::Float64
    sum = 0
    dim = int(sqrt(length(mat)))
    for i in length(mat)
        dn = i / dim
        ds = i % dim
        sum += mat[i]*G.pdf(dN, dS)
    end
    return sum
end

function logLikelihood(A::Array{Int8}, mat::Array{Float64}, G::Distribution)::Float64
    sum = 0
    for Ai in A
        sum += log(likelihood(Ai, mat, G))
    end
end

function discretize(G::Distribution, inputSpace::Array{Float64})::Distribution
    dict = {}
    normalization = 0
    for input in inputSpace
        pdf = G.pdf(input)
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