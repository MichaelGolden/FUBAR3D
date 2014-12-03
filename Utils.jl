function copy{T}(source::Array{T, 1}, dest::Array{T, 1})
    n = length(source)
    for i = 1:n
        dest[i] = source[i]
    end
end

function copy{T}(source::Array{T, 2}, dest::Array{T, 2})
    rows = size(source, 1)
    cols = size(source, 2)
    for col=1:cols, row=1:rows
        dest[row, col] = source[row, col]
    end
end

function logSumExp(a::Float64, b::Float64)
    m = max(a,b)
    return m + log(exp(a-m)+exp(b-m))
end

function logSumOneAndExponential(a::Float64)
    if a < 0.0
        return log1p(exp(a))
    else
        return a + log1p(exp(-a))
    end
end

function sample(rng, v::Array{Float64})
    s = sum(v)
    n = length(v)
    r = rand(rng)
    cumsum = 0.0
    for i=1:n
        cumsum += v[i] / s
        if cumsum >= r
            return i
        end
    end

    return n
end

nonzeros(S::SparseMatrixCSC) = S.nzval
rowvals(S::SparseMatrixCSC) = S.rowval
nzrange(S::SparseMatrixCSC, col::Integer) = S.colptr[col]:(S.colptr[col+1]-1)

