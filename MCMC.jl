using Distributions

import Base.copy
include("Utils.jl")
include("MCMCState.jl")
include("AcceptanceLogger.jl")

function computeLogLikelihoods(state::MCMCState, conditionals::Array{Float64,2}, scalarSum::Float64, c::Float64, partialPhyloUpdate::Bool, partialMRFUpdate::Bool)
    if !partialPhyloUpdate # recompute site likelihoods
        state.cachedPhyloLikelihoods = log((state.θ)'*conditionals)
    end

    phyloLikelihood = 0.0
    for i = 1:numSites
        phyloLikelihood += state.cachedPhyloLikelihoods[state.hiddenStates[i],i]
    end

    state.logPhyloLikelihood = phyloLikelihood+scalarSum
    state.logMRFLikelihood = computeAutoLogisticMRFLogLikelihood(state)
    state.logMRFLikelihood = computeGaussianMRFLogLikelihood(state, partialMRFUpdate)
    state.logLikelihood = state.logPhyloLikelihood+state.logMRFLikelihood

    n=size(state.cachedPhyloLikelihoods, 1) # get number of hidden states
    state.logPrior = 0.0
    for hiddenState=1:n
        state.logPrior += sum(log(state.θ[:,hiddenState])*(c-1.0))
    end
end

function gibbsUpdate(state::MCMCState, hiddenState::Int, conditionals::Array{Float64,2}, scalarSum::Float64, c::Float64)
    φ = zeros(Float64, numGridPoints)
    v = zeros(Float64, numGridPoints)
    for i=1:numSites
        if state.hiddenStates[i] == hiddenState
            for j=1:numGridPoints
                v[j] = state.θ[j,hiddenState]*conditionals[j,i]
            end
            φ[sample(v)] += 1.0
        end
    end
    alphas = zeros(Float64, numGridPoints)
    for j=1:numGridPoints
        alphas[j] = φ[j]+c
    end

    d = Dirichlet(alphas)
    r = rand(d)
    for i=1:numGridPoints
        state.θ[i,hiddenState] = r[i]
    end
end

function fractionNeighboursSame(hiddenStates::Array{Int8})
    n = length(hiddenStates)
    fractionSame = 0.0
    for i=1:n
        same = 0
        total = 0
        for j=i-1:i+1
            if i != j && j >= 1 && j <= n
                if hiddenStates[i] == hiddenStates[j]
                    same += 1
                end
                total += 1
            end
        end
        fractionSame += same / total
    end

    return fractionSame / n
end

function thetaMove(θ::Array{Float64})
    x = rand(1:numGridPoints)
    y = x
    while x==y
        y = rand(1:numGridPoints)
    end
    η = (1.0/numGridPoints)*rand(rng)*2.0
    θ[x] += η
    θ[y] -= η

    if θ[y] < 0.0 || θ[x] > 1.0
        return false
    else
        θ = (1.0/sum(θ))*θ # normalise θ to sum to 1
        return true
    end
end

function computeAutoLogisticMRFLogLikelihood(state::MCMCState)
#= computes auto-logistic likelihood where neighbours are -3 to +3 adjacent positions in sequence =#
    ll = 0.0
    dist = 3
    for i=1:numSites # sites
        x_i = state.hiddenStates[i]-1.0
        denominator = 0.0
        start = max(1,i-dist)
        stop = min(numSites,i+dist)
        for j=start:stop # j = neighbours of i
            x_j = state.hiddenStates[j]-1.0
            if i != j
                mult = state.β*exp(-abs(i-j-1)*state.τ) # i-j-1 starts at 0
                if x_i == 0.0
                    denominator += state.α*x_i + mult*(2.0*x_j-1.0)
                else
                    denominator -= state.α*x_i + mult*(2.0*x_j-1.0)
                end
            end
        end
        cond = -logSumOneAndExponential(denominator)
        ll += cond
    end
    return ll
end

function computeGaussianMRFLogLikelihood(state::MCMCState, partialMRFUpdate::Bool)
#= computes gaussian likelihood where neighbours are -3 to +3 adjacent positions in sequence =#
    B=zeros(Float64, numSites, numSites)

    dist = 3
    for i=1:numSites
        start = max(1,i-dist)
        stop = min(numSites,i+dist)
        for j=start:stop # j = neighbours of i
            if i == j
                B[i,j] = 1.0
            else
                B[i,j] = state.β*exp(-abs(i-j-1)*state.τ) # i-j-1 starts at 0
                B[j,i] = B[i,j]
            end
        end
    end

    if !partialMRFUpdate # recompute determinant: O(N^3)
        state.cachedDet = logdet(B)
    end

    sigma = 1.0 # note I've fixed sigma, this should probably be a free parameter
    v = zeros(Float64, numSites)
    for i=1:numSites
        v[i] = state.hiddenStates[i] - 1.5
    end

    ll = (state.cachedDet*0.5) - (log(2*pi*sigma*sigma)*numSites*0.5) + (((v')*B*v)[1]/(2.0*sigma*sigma))

    return ll
end

gridInfoFile = "datasets/hiv1_env_300.nex.grid_info"
#gridInfoFile = "datasets/hcv1_polyprotein_300.nex.grid_info"
#gridInfoFile = "datasets/lysin.nex.grid_info"

srand(948402288028201)
rng = MersenneTwister(948402288028201)
c=0.5

include("FUBARDataset.jl")
dataset = loadDataset(gridInfoFile)
numGridPoints = dataset[1]
numSites = dataset[2]
grid = dataset[3]
conditionals = dataset[4]
scalars = dataset[5]
scalarSum = sum(scalars)

numHiddenStates=2
partialPhyloUpdate = false
partialMRFUpdate = false

initialθ = zeros(Float64, numGridPoints, numHiddenStates)
for col=1:numHiddenStates
    initialθ[:,col] = rand(numGridPoints)
    initialθ[:,col] = initialθ[:,col]/sum(initialθ[:,col])
end
initialPhyloLikelihoods = zeros(Float64, numHiddenStates, numSites)

currentState = MCMCState(0.0,0.0,0.0,0.0,copy(initialθ),[rand(1:2) for i=1:numSites],0.0,0.5,1.0,copy(initialPhyloLikelihoods),0.0)
computeLogLikelihoods(currentState, conditionals, scalarSum, c, partialPhyloUpdate, partialMRFUpdate)
proposedState  = MCMCState(currentState)

iterations = 0
logWriter = open("mcmc.log","w")
thetaWriter = open("theta.log","w")
thetaWriter2 = open("theta2.log","w")
hiddenWriter = open("hidden.log","w")
acceptanceWriter = open("acceptance.log","w")

moveWeights = [1.0, 250.0, 25.0]
#moveWeights = [1.0, 0.0]

logger = AcceptanceLogger()
moveDescription = ""
firstIteration = true

sampleRate=100

maxIterations=100000000
for iter=1:maxIterations
    valid = true
    gibbs = false
    move = sample(moveWeights)

    if move == 1
        #valid = valid ? thetaMove(proposedState.θ[:,rand(1:2)] : false
        gibbsUpdate(proposedState, 1, conditionals, scalarSum, c)
        gibbsUpdate(proposedState, 2, conditionals, scalarSum, c)
        gibbs=true
        partialPhyloUpdate = false
        partialMRFUpdate = true
        moveDescription = "gibbs"
    elseif move == 2 # hidden state move
        n = rand(1:5)
        for i=1:n
            proposedState.hiddenStates[rand(1:numSites)] = rand(1:2)
        end
        partialPhyloUpdate = true
        partialMRFUpdate = true
        moveDescription = "hiddenStates_$n"
    elseif move == 3 # parameter moves
        proposedState.α += randn(rng) * 0.1
        proposedState.β =  proposedState.β*exp(randn(rng) * 0.2)
        proposedState.τ += randn(rng) * 0.1
        if proposedState.β <= 0 || proposedState.τ <= 0 || abs(proposedState.β) >= 5
            valid = false
        end
        partialPhyloUpdate = false
        partialMRFUpdate = false
        moveDescription = "parameters"
    end
    proposedState.α = 0.0
    #proposedState.τ = 0.0
    #proposedState.β = 0.0

    logPropRatio = 0
    if valid        
        computeLogLikelihoods(proposedState, conditionals, scalarSum, c, partialPhyloUpdate, partialMRFUpdate)
        logPropRatio =  proposedState.logPrior - currentState.logPrior
    end
    Δll = proposedState.logLikelihood - currentState.logLikelihood + logPropRatio
    if valid && (gibbs || exp(Δll) > rand(rng))
        copy(proposedState, currentState)
        logAccept!(logger, moveDescription)
    else
        copy(currentState, proposedState)
        logReject!(logger, moveDescription)
    end

    if iter % sampleRate == 0
        if firstIteration
            write(logWriter, "iter\tLL\tphyloLL\tMRFLL\tpriorLL\tfracNeighboursSame\talpha\tbeta\ttau\n")
            firstIteration = false
        end
        write(logWriter,string(iter, "\t", currentState.logLikelihood,"\t", currentState.logPhyloLikelihood,"\t", currentState.logMRFLikelihood,"\t", currentState.logPrior, "\t", fractionNeighboursSame(currentState.hiddenStates), "\t", currentState.α, "\t", currentState.β, "\t", currentState.τ, "\n"))
        flush(logWriter)

        if iter % (sampleRate*10) == 0
            write(thetaWriter,join(currentState.θ[:,2],"\t"),"\n")
            write(thetaWriter2,join(currentState.θ[:,2],"\t"),"\n")
            flush(thetaWriter)
            flush(thetaWriter2)
            write(hiddenWriter,join(currentState.hiddenStates,""), "\n")
            flush(hiddenWriter)
        end

        if iter % (sampleRate*100) == 0
            write(acceptanceWriter,iter, join(list(logger),"\t"),"\n")
            flush(acceptanceWriter)
            clear!(logger)
        end
    end
end
close(logWriter)
close(thetaWriter)
close(thetaWriter2)
close(hiddenWriter)
close(acceptanceWriter)
