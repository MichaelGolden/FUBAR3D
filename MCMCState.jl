include("Utils.jl")

type MCMCState
    logLikelihood::Float64
    logPhyloLikelihood::Float64
    logMRFLikelihood::Float64
    logPrior::Float64
    θ::Array{Float64, 2}
    hiddenStates::Array{Int8}
    α::Float64
    β::Float64
    τ::Float64
    cachedPhyloLikelihoods::Array{Float64,2}

    MCMCState(logLikelihood,logPhyloLikelihood,logMRFLikelihood,logPrior,θ,hiddenStates,α,β,τ,cachedPhyloLikelihoods) = new(logLikelihood,logPhyloLikelihood,logMRFLikelihood,logPrior,θ,hiddenStates,α,β,τ,cachedPhyloLikelihoods)

    function MCMCState(other::MCMCState)
        return new(other.logLikelihood,other.logPhyloLikelihood,other.logMRFLikelihood,other.logPrior,copy(other.θ),copy(other.hiddenStates),other.α,other.β,other.τ,copy(other.cachedPhyloLikelihoods))
    end
end

function copy(source::MCMCState, dest::MCMCState)
    dest.logLikelihood = source.logLikelihood
    dest.logPhyloLikelihood = source.logPhyloLikelihood
    dest.logMRFLikelihood = source.logMRFLikelihood
    dest.logPrior = source.logPrior
    copy(source.θ, dest.θ)
    copy(source.hiddenStates, dest.hiddenStates)
    dest.α = source.α
    dest.β = source.β
    dest.τ = source.τ
    copy(source.cachedPhyloLikelihoods, dest.cachedPhyloLikelihoods)
end
