
######
# Define multivariate distribution
######

using Distributions
using Plotly


## define multivariate distribution with rho
  immutable MultiVnorm <: AbstractMvNormal
      μ::Array{Float64,1}
      σ::Array{Float64,1}
      ρ::Float64
      #constructor
      MultiVnorm(μ,σ,ρ) = new(μ,σ,ρ)
  end

  # calculate pdf on mvn
  function pdf(d::MultiVnorm,x::Array{Float64,1})
    Σ = [d.σ[1]^2 d.ρ*d.σ[1]*d.σ[2] ; d.ρ*d.σ[1]*d.σ[2] d.σ[2]]
    C = 1/sqrt(((2*π)^2)*det(Σ))
    E = exp(-0.5*((x-d.μ)')*(Σ^-1)*(x-d.μ))
    (C*E)[1]
#    println(C*E)
#    p = float64(C*E)[1]
#    println(p)
#    println(typeof(p))
 #   p
end

function test_BiVarNorm()
  x=1:.01:5
  y=1:.01:5
  d = MultiVnorm([2.5,2.5],[1.,1.],.1)
  #pdf = [BiVnormPDF([i,j],[2.5,2.5],[1.,1.],.1) for i=x, j=y]
  P = [pdf(d,[i,j]) for i=x, j=y]
  visualize(P,x,y)
end

function visualize(pdf,x,y)
  Plotly.signin("bkellman", "49ez4efxbk")
  data = [[
      "z" => pdf,
      "x" => x,
      "y" => y,
      "type" => "heatmap"]]
  response = Plotly.plot(data, ["BiVariateNormal" => "basic-heatmap", "fileopt" => "overwrite"])
  println(response["url"])
end

test_BiVarNorm()


#function BiVnormPDF(x::Array{Float64,1},μ::Array{Float64,1},σ::Array{Float64,1},ρ::Float64)
#  Σ = [σ[1]^2 ρ*σ[1]*σ[2] ; ρ*σ[1]*σ[2] σ[2]]
#  BiVnormPDF(x,μ,Σ)
#end

#function BiVnormPDF(x::Array{Float64,1},μ::Array{Float64,1},Σ::Array{Float64,2})
#  C = 1/sqrt(((2*π)^2)*det(Σ))
#  E = exp(-0.5*((x-μ)')*(Σ^-1)*(x-μ))
#  C*E
#end

#  # calculate pdf on mvn
#  function pdf(d::MultiVnorm,x,y)
#    coef = 1/(2*π*d.σ1*d.σ2*sqrt(1-ρ^2))
#
#  end
#end

#function BiVnormPDF(x::Array{Float64,1},μ1::Float64,μ2::Float64,σ1::Float64,σ2::Float64,ρ::Float64)
#  μ = [μ1,μ2]
#  Σ = [σ1^2 ρ*σ1*σ2 ; ρ*σ1*σ2 σ2]
#  BiVnormPDF(x,μ,Σ)
#end
