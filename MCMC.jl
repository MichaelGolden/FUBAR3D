using Distributions

import Base.copy
include("Utils.jl")
include("MCMCState.jl")
include("AcceptanceLogger.jl")

function computeLogLikelihoods(numGridPoints::Int, numSites::Int, state::MCMCState, conditionals::Array{Float64,2}, scalarSum::Float64, c::Float64, distMatrix::SparseMatrixCSC{Float64, Int}, B::Array{Float64,2}, partialPhyloUpdate::Bool, partialMRFUpdate::Bool)
    if !partialPhyloUpdate # recompute site likelihoods
        state.cachedPhyloLikelihoods = log((state.θ)'*conditionals)
    end

    phyloLikelihood = 0.0
    for i = 1:numSites
        phyloLikelihood += state.cachedPhyloLikelihoods[state.hiddenStates[i],i]
    end

    state.logPhyloLikelihood = phyloLikelihood+scalarSum
    state.logPhyloLikelihood = 0.0 # no data
    state.logMRFLikelihood = computeIsingPseudoLogLikelihood(numGridPoints, numSites, state, distMatrix)
    #state.logMRFLikelihood = computeBensMRFLogLikelihood(numGridPoints, numSites, state, distMatrix)
    #state.logMRFLikelihood = computeGaussianMRFLogLikelihood(numGridPoints, numSites, state, distMatrix, B, partialMRFUpdate)

    state.logLikelihood = state.logPhyloLikelihood+state.logMRFLikelihood

    n=size(state.cachedPhyloLikelihoods, 1) # get number of hidden states
    state.logPrior = 0.0
    for hiddenState=1:n
        state.logPrior += sum(log(state.θ[:,hiddenState])*(c-1.0))
    end
    state.logPrior = 0.0  # no data
end

function gibbsUpdate(rng::MersenneTwister,numGridPoints::Int, numSites::Int, state::MCMCState, hiddenState::Int, conditionals::Array{Float64,2}, scalarSum::Float64, c::Float64)
    φ = zeros(Float64, numGridPoints)
    v = zeros(Float64, numGridPoints)
    for i=1:numSites
        if state.hiddenStates[i] == hiddenState
            for j=1:numGridPoints
                v[j] = state.θ[j,hiddenState]*conditionals[j,i]
            end
            φ[sample(rng, v)] += 1.0
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

function thetaMove(numGridPoints::Int, numSites::Int, θ::Array{Float64}) # unused
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

function computeIsingPseudoLogLikelihood(numGridPoints::Int, numSites::Int, state::MCMCState, distMatrix::SparseMatrixCSC{Float64, Int})
    ll = 0.0
    rows = rowvals(distMatrix)
    vals = nonzeros(distMatrix)
    for i=1:numSites
        x_i = state.hiddenStates[i]-1.0
        q = 0.0
        for sparseindex in nzrange(distMatrix, i)
            j = rows[sparseindex]  # j = neighbours of i
            mult = state.β*exp(-(vals[sparseindex]-1.0)*state.τ) # i-j-1 starts at 0
            x_j = state.hiddenStates[j]-1.0
            if x_i == x_j
                q += mult
            else
                q -= mult
            end
        end
        ll += -logSumOneAndExponential(q)
    end
    return ll
end

function computeBensMRFLogLikelihood(numGridPoints::Int, numSites::Int, state::MCMCState, distMatrix::SparseMatrixCSC{Float64, Int})
    ll = 0.0
    rows = rowvals(distMatrix)
    vals = nonzeros(distMatrix)
    for i=1:numSites # sites
        x_i = state.hiddenStates[i]-1.0
        v = 0.0
        k = 0.0
        for sparseindex in nzrange(distMatrix, i)
            j = rows[sparseindex]
            v += state.hiddenStates[j]-1.0
            k += 1.0
        end
        t_eta = (1.0/(1.0-state.β))-1.0
        #print(state.β, "\t", t_eta, "\n")
        numerator = t_eta + v
        denominator = k + 2.0*t_eta
        cond = log(numerator)-log(denominator)
        ll += cond
    end
    return ll
end


function computeGaussianMRFLogLikelihood(numGridPoints::Int, numSites::Int, state::MCMCState, distMatrix::SparseMatrixCSC{Float64, Int}, B::Array{Float64,2}, partialMRFUpdate::Bool)
    B=zeros(Float64, numSites, numSites)

    rows = rowvals(distMatrix)
    vals = nonzeros(distMatrix)
    for i=1:numSites
        B[i,i] = 1.0 * (state.σ*state.σ)
        for sparseindex in nzrange(distMatrix, i) # j = neighbours of i
            j = rows[sparseindex]
            B[i,j] = state.β * (state.σ*state.σ) * exp(-(vals[sparseindex]-1)*state.τ)
            B[j,i] = B[i,j]
        end
    end

    B=B'B

    #writer=open("matrix.log","w")
   # write(writer,string(B),"\n")
    #close(writer)

    if !isposdef(B)
        print("Not positive definite", state.β,"\n")
        return -Inf
    end

    if !partialMRFUpdate # recompute determinant: O(N^3)
        state.cachedDet = logdet(B)
    end

    v = zeros(Float64, numSites)

    for i=1:numSites
        v[i] = state.hiddenStates[i] - 1.5
    end


    ll = -(2*pi*numSites*0.5) - (state.cachedDet*0.5) - (0.5*(v'*inv(B)*v)[1])

    return ll
end

function computeGaussianMRFLogLikelihoodOld(numGridPoints::Int, numSites::Int, state::MCMCState, distMatrix::SparseMatrixCSC{Float64, Int}, partialMRFUpdate::Bool)
    B=zeros(Float64, numSites, numSites)

    rows = rowvals(distMatrix)
    vals = nonzeros(distMatrix)
    for i=1:numSites
        B[i,i] = 1.0
        for sparseindex in nzrange(distMatrix, i) # j = neighbours of i
            j = rows[sparseindex]
            B[i,j] = state.β*exp(-(vals[sparseindex]-1)*state.τ)
            B[j,i] = B[i,j]
        end
    end

    if !isposdef(B)
        return -Inf
    end

    if !partialMRFUpdate # recompute determinant: O(N^3)
        state.cachedDet = logdet(B)
    end

    v = state.hiddenStates[i] - 1.5 # element-wise subtraction

    a = (state.cachedDet*0.5)
    b = (log(2*pi*state.σ*state.σ)*numSites*0.5)
    c = (((v')*B*v)[1]/(2.0*state.σ*state.σ))
    print(v,"\t",state.σ,"\t",state.cachedDet,"\t",((v')*B*v),"\t",a,"\t",b,"\t",c,"\n")

    ll = (state.cachedDet*0.5) - (log(2*pi*state.σ*state.σ)*numSites*0.5) + (((v')*B*v)[1]/(2.0*state.σ*state.σ))

    return ll
end

function getLinearDependenceDistanceMatrix(numSites::Int, dist::Int)
    I = Int[]
    J = Int[]
    V = Float64[]

    for i=1:numSites
        start = max(1,i-dist)
        stop = min(numSites,i+dist)
        for j=start:stop # j = neighbours of i
            if i != j
                push!(I, i)
                push!(J, j)
                push!(V, abs(i-j))
            end
        end
    end
    return sparse(I, J, V)
end

function mainMCMC(distMatrix::SparseMatrixCSC{Float64,Int64} , gridInfoFile , β::Float64=0 ,maxIterations=1e4,repNum=0)
    #gridInfoFile = "datasets/hiv1_env_300.nex.grid_info"
    #gridInfoFile = "datasets/hcv1_polyprotein_300.nex.grid_info"
    #gridInfoFile = "datasets/lysin.nex.grid_info"

    srand(948402288028201+repNum)
    rng = MersenneTwister(948402288028201+repNum)
    c=0.5

    include("FUBARDataset.jl")
    dataset = loadDataset(gridInfoFile)
    numGridPoints = dataset[1]
    numSites = dataset[2]
    grid = dataset[3]
    conditionals = dataset[4]
    scalars = dataset[5]
    scalarSum = sum(scalars)

    #distMatrix = getLinearDependenceDistanceMatrix(numSites, 1)

    #distMatrix = sparse(readdlm("distance.matrix.thresh15.csv",' '))

    numHiddenStates=2
    partialPhyloUpdate = false
    partialMRFUpdate = false

    initialθ = zeros(Float64, numGridPoints, numHiddenStates)
    for col=1:numHiddenStates
        initialθ[:,col] = rand(numGridPoints)
        initialθ[:,col] = initialθ[:,col]/sum(initialθ[:,col])
    end
    initialPhyloLikelihoods = zeros(Float64, numHiddenStates, numSites)

    B=zeros(Float64, numSites, numSites)

    currentState = MCMCState(0.0,0.0,0.0,0.0,copy(initialθ),[rand(1:2) for i=1:numSites],0.0,0.5,1.0,1.0,copy(initialPhyloLikelihoods),0.0)
    computeLogLikelihoods(numGridPoints, numSites, currentState, conditionals, scalarSum, c, distMatrix, B, partialPhyloUpdate, partialMRFUpdate)
    proposedState  = MCMCState(currentState)

    iterations = 0

    fnamesary = ["mcmc.log","theta.log","theta2.log","hidden.log","acceptance.log"]
    gridInfoLogFile = split(gridInfoFile,"/")[2]
    fnamesary = ["log/"*gridInfoLogFile*"_"*string(fnamesary[i]) for i=1:length(fnamesary)]
#    fnamesary = [string(gridInfoLogFile,"_",string(fnamesary[i])) for i=1:length(fnamesary)]
    fi=0
    logWriter = open(fnamesary[fi+=1],"w")
    thetaWriter = open(fnamesary[fi+=1],"w")
    thetaWriter2 = open(fnamesary[fi+=1],"w")
    hiddenWriter = open(fnamesary[fi+=1],"w")
    acceptanceWriter = open(fnamesary[fi+=1],"w")

    #logWriter = open("mcmc.log","w")
    #thetaWriter = open("theta.log","w")
    #thetaWriter2 = open("theta2.log","w")
    #hiddenWriter = open("hidden.log","w")
    #acceptanceWriter = open("acceptance.log","w")

    if β > 0
      moveWeights = [1.0, 250.0, 0]  # set move 3 from 20 -> 0
    else
      moveWeights = [1.0, 250.0, 20]  # set move 3 from 20 -> 0
    end
    #moveWeights = [1.0, 0.0]

    logger = AcceptanceLogger()
    moveDescription = ""
    firstIteration = true


    sampleRate=100

    #maxIterations=100000
    for iter=1:maxIterations
        valid = true
        gibbs = false
        move = sample(rng, moveWeights)

        if move == 1
            #valid = valid ? thetaMove(proposedState.θ[:,rand(1:2)] : false
            gibbsUpdate(rng, numGridPoints, numSites, proposedState, 1, conditionals, scalarSum, c)
            gibbsUpdate(rng, numGridPoints, numSites, proposedState, 2, conditionals, scalarSum, c)
            gibbs=true
            partialPhyloUpdate = false
            partialMRFUpdate = true
            moveDescription = "gibbs"
        elseif move == 2 # hidden state move
            n = rand(1:4)
            for i=1:n
                proposedState.hiddenStates[rand(1:numSites)] = rand(1:2)
            end
            partialPhyloUpdate = true
            partialMRFUpdate = true
            moveDescription = "hiddenStates_$n"
        elseif move == 3 # parameter moves
            proposedState.α += randn(rng) * 0.1
            proposedState.β +=  randn(rng) * 0.25
            #proposedState.β +=  randn(rng) * 0.02
            proposedState.σ +=  randn(rng) * 0.01
            proposedState.τ +=  randn(rng) * 0.1
            if proposedState.β < 0.0  || proposedState.β >= 5.0 || proposedState.σ  <= 0.0|| proposedState.σ > 10 || proposedState.τ < 0.0
                valid = false
            end
            partialPhyloUpdate = false
            partialMRFUpdate = false
            moveDescription = "parameters"
        end
        proposedState.α = 0.0
        proposedState.τ = 0.0
        #proposedState.σ = 1.0
        proposedState.β = β  # fix beta here

        logPropRatio = 0
        if valid
            computeLogLikelihoods(numGridPoints, numSites, proposedState, conditionals, scalarSum, c, distMatrix, B, partialPhyloUpdate, partialMRFUpdate)
            #prior = log(proposedState.β)*3.0
            prior = 0.0
            logPropRatio =  proposedState.logPrior - currentState.logPrior + prior
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
                write(logWriter, "iter\tLL\tphyloLL\tMRFLL\tpriorLL\tfracNeighboursSame\talpha\tbeta\tsigma\ttau\n")
                firstIteration = false
            end
            write(logWriter,string(iter, "\t", currentState.logLikelihood,"\t", currentState.logPhyloLikelihood,"\t", currentState.logMRFLikelihood,"\t", currentState.logPrior, "\t", fractionNeighboursSame(currentState.hiddenStates), "\t", currentState.α, "\t", currentState.β,"\t", currentState.σ, "\t", currentState.τ, "\n"))
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
                write(acceptanceWriter,string(iter),"\t", join(list(logger),"\t"),"\n")
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

    return currentState.hiddenStates
end

#main()

