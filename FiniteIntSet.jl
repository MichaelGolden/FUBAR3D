type FiniteIntSet
    maxElements::Int
    array::Array{Int}
    contains::Array{Int}
    size::Int
    state::Int

    MCMCState(maxElements, array, contains, size, state) = new(maxElements, array, contains, size, state)

    function FiniteIntSet(maxElements::Int)
        array = zeros(Int,maxElements)
        contains = zeros(Int,maxElements)
        size = 0
        state = 1
        return new(maxElements, array, contains, size, state)
    end
end

function push!(set::FiniteIntSet, val::Int)
    if(set.contains[val] != set.state)
        set.array[set.size+1] = val
        set.size += 1
    end
    set.contains[val] = set.state
end

function clear!(set::FiniteIntSet)
    set.state += 1
    set.size = 0
end

function contains(set::FiniteIntSet, val::Int)
    return set.contains[val] == set.state
end

intset = FiniteIntSet(10)
push!(intset, 8)
print(contains(intset,5),"\n")
print(contains(intset,8),"\n")
push!(intset, 2)
clear!(intset)
push!(intset, 5)
print(contains(intset,5),"\n")
print(contains(intset,8),"\n")
