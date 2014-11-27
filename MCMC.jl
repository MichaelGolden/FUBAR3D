using Distributions

import Base.copy
include("Utils.jl")
include("MCMCState.jl")

#gridInfoFile = "datasets/hiv1_env_300.nex.grid_info"
#gridInfoFile = "datasets/hcv1_polyprotein_300.nex.grid_info"
gridInfoFile = "datasets/lysin.nex.grid_info"

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

function computeLogLikelihoods(state::MCMCState, conditionals::Array{Float64,2}, scalarSum::Float64, c::Float64)
    siteLikelihood = log((state.θ)'*conditionals)
    phyloLikelihood = 0.0
    for i = 1:numSites
        phyloLikelihood += siteLikelihood[state.hiddenStates[i],i]
    end

    state.logPhyloLikelihood = phyloLikelihood+scalarSum
    state.logMRFLikelihood = computeAutobinomialMRFLogLikelihood(state.hiddenStates, state.α, state.β,  state.τ)
    state.logLikelihood = state.logPhyloLikelihood+state.logMRFLikelihood

    n=size(siteLikelihood, 1) # get number of hidden states
    state.logPrior = 0.0
    for hiddenState=1:n
        state.logPrior += sum(log(state.θ[:,hiddenState])*(c-1.0))
    end
end

function gibbsStep(state::MCMCState, hiddenState::Int, conditionals::Array{Float64,2}, scalarSum::Float64, c::Float64)
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

function fractionSameNeighbours(hiddenStates::Array{Int8})
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

function computeAutobinomialMRFLogLikelihood(hiddenStates::Array{Int8}, α::Float64, β::Float64, τ::Float64)
    ll = 0.0
    dist = 3
    for i=1:numSites
        x_i = hiddenStates[i]-1
        numerator = α*x_i
        denominator = α
        start = max(1,i-dist)
        stop = min(numSites,i+dist)
        for j=start:stop
            x_j = hiddenStates[j]-1
            if i != j
                mult = β*exp(-abs(i-j-1)*τ) # i-j-1 starts at 0
                numerator += mult*x_i*x_j
                denominator += mult*x_j
            end
        end
        cond = numerator - logSumOneAndExponential(denominator)
        ll += cond
    end
    return ll
end

function computeAutobinomialMRFLogLikelihood2(hiddenStates::Array{Int8}, α::Float64, β::Float64, τ::Float64)
    ll = 0.0
    dist = 3
    for i=1:numSites
        x_i = hiddenStates[i]-1.0
        denominator = 0.0
        start = max(1,i-dist)
        stop = min(numSites,i+dist)
        for j=start:stop
            x_j = hiddenStates[j]-1.0
            if i != j
                mult = β*exp(-abs(i-j-1)*τ) # i-j-1 starts at 0
                #numerator += mult*x_i*x_j
                if x_i == 0.0
                    denominator += α*x_i + mult*(2.0*x_j-1.0)
                else
                    denominator -= α*x_i + mult*(2.0*x_j-1.0)
                end
            end
        end
        if i == 1
            print(denominator,"\n")
        end
        cond = -logSumOneAndExponential(denominator)
        ll += cond
    end
    return ll
end

function computeGaussianMRFLogLikelihood(hiddenStates::Array{Int8}, α::Float64, β::Float64, τ::Float64)
    B=zeros(Float64, numSites, numSites)

    for i=1:numSites
        for j=1:numSites
            if i != j
                B[i,j] = β*exp(-abs(i-j-1)*τ) # i-j-1 starts at 0
            else
                B[i,j] = 1.0
            end
        end
    end

    logd = logdet(B)
    sigma = 1.0
    v = zeros(Float64, numSites)
    for i=1:numSites
        v[i] = hiddenStates[i] - 1.5
    end

    ll = (logd*0.5) - (log(2*pi*sigma*sigma)*numSites*0.5) + (((v')*B*v)[1]/(2.0*sigma*sigma))

    return ll
end

numHiddenStates=2
initialθ = zeros(Float64, numGridPoints, numHiddenStates)
for col=1:numHiddenStates
    initialθ[:,col] = rand(numGridPoints)
    initialθ[:,col] = initialθ[:,col]/sum(initialθ[:,col])
end

currentState = MCMCState(0.0,0.0,0.0,0.0,copy(initialθ),[rand(1:2) for i=1:numSites],0.0,0.5,1.0)
computeLogLikelihoods(currentState, conditionals, scalarSum, c)
proposedState  = MCMCState(currentState)

iterations = 0
logWriter = open("mcmc.log","w")
thetaWriter = open("theta.log","w")
thetaWriter2 = open("theta2.log","w")
hiddenWriter = open("hidden.log","w")
weights = [1.0, 250.0, 25.0]
#weights = [1.0, 0.0]
include("AcceptanceLogger.jl")
logger = AcceptanceLogger()
moveDescription = ""
maxIterations=100000000
for iter=1:maxIterations
    valid = true
    gibbs = false
    move = sample(weights)
    if move == 1
        #valid = valid ? thetaMove(proposedState.θ[:,rand(1:2)] : false
        gibbsStep(proposedState, 1, conditionals, scalarSum, c)
        gibbsStep(proposedState, 2, conditionals, scalarSum, c)
        gibbs=true
        moveDescription = "gibbs"
    elseif move == 2
        n = rand(1:5)
        for i=1:n
            proposedState.hiddenStates[rand(1:numSites)] = rand(1:2)
        end
        moveDescription = "hiddenStates_$n"
    elseif move == 3
        proposedState.α += randn(rng) * 0.1
        proposedState.β =  proposedState.β*exp(randn(rng) * 0.2)
        proposedState.τ += randn(rng) * 0.1
        if proposedState.β <= 0 || proposedState.τ <= 0 || abs(proposedState.β) >= 5
            valid = false
        end
        moveDescription = "parameters"
    end
    proposedState.α = 0.0
    #proposedState.τ = 0.0
    #proposedState.β = 10000.0
    #=
    proposedState.hiddenStates = [1 for i=1:numSites]
    mrfll = computeAutobinomialMRFLogLikelihood2(proposedState.hiddenStates, proposedState.α, proposedState.β,  proposedState.τ)
    mrfll2 = computeGaussianMRFLogLikelihood(proposedState.hiddenStates, proposedState.α, proposedState.β,  proposedState.τ)
    print("A",mrfll,"\t",mrfll2,"\t", proposedState.β,"\n")
    proposedState.hiddenStates = [2 for i=1:numSites]
    mrfll = computeAutobinomialMRFLogLikelihood2(proposedState.hiddenStates, proposedState.α, proposedState.β,  proposedState.τ)
    mrfll2 = computeGaussianMRFLogLikelihood(proposedState.hiddenStates, proposedState.α, proposedState.β,  proposedState.τ)
    print("B",mrfll,"\t",mrfll2,"\t", proposedState.β,"\n")
    =#

    logPropRatio = 0
    if valid        
        computeLogLikelihoods(proposedState, conditionals, scalarSum, c)
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

    if iter % 50 == 1
        if iter == 1
            write(logWriter, "iter\tLL\tphyloLL\tMRFLL\tpriorLL\tfractionSame\talpha\tbeta\ttau\n")
        end
        fractionSame = fractionSameNeighbours(currentState.hiddenStates)
        write(logWriter,string(iter, "\t", currentState.logLikelihood,"\t", currentState.logPhyloLikelihood,"\t", currentState.logMRFLikelihood,"\t", currentState.logPrior, "\t", fractionSame, "\t", currentState.α, "\t", currentState.β, "\t", currentState.τ, "\n"))
        flush(logWriter)
        write(thetaWriter,string(currentState.θ[:,1], "\n"))
        write(thetaWriter2,string(currentState.θ[:,2], "\n"))
        flush(thetaWriter)
        flush(thetaWriter2)
        write(hiddenWriter,string(currentState.hiddenStates, "\n"))
        flush(hiddenWriter)

        if iter % 5000 == 1
            print(list(logger),"\n")
            clear!(logger)
        end
        #=
        if iterations % 500 == 0
            gll = computeGaussianMRFLogLikelihood(currentState.hiddenStates, currentState.α, currentState.β,  currentState.τ)
            print(currentState.logMRFLikelihood, "\t", gll, "\n")
        end
        =#
    end
end
close(logWriter)
close(thetaWriter)
close(thetaWriter2)
close(hiddenWriter)
