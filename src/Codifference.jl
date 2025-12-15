module Codifference

using StatsBase, Distributions

export ecf, lcf, cod, lcfAsymptDistr, codAsymptDistr, lcfConfInterval, codConfInterval


"""
    ecf(X::AbstractVector, θ::Real)

Empirical characteristic function for random vectors with symmetric distribution. Returns real valued output
"""
ecf(X::AbstractVector, θ::Real) = mean( cos(θ*x) for x in X)

"""
    lcf(X::AbstractVector, θ::Real = 1)

Empirical log characteristic function of sample X. Returns NaN if emprical characteristic function is non-positive.
"""
function lcf(X::AbstractVector, θ::Real = 1) 
    θ ≈ 0 && return one(θ)*var(X) 
    e = ecf(X,θ)
    e <= 0 && return one(θ)*eltype(X)(NaN)
    return -2/θ^2 * log(e)
end

"""
    cod(X::AbstractVector, Y::AbstractVector, θ::Real = 1, type::Symbol = :s)

Empirical codifference of sample X, Y.
"""
function cod(X::AbstractVector, Y::AbstractVector, θ::Real = 1, type::Symbol = :s)
    type ∉ (:s, :+, :-) && error("Unrecognised type of codifference")
    if type == :s
        1/4 * ( lcf(X+Y, θ) - lcf(X-Y, θ) )
    elseif type == :+
        1/2 * ( lcf(X, θ) + lcf(Y, θ) - lcf(X-Y, θ) )
    else
        -1/2 * ( lcf(X, θ) + lcf(Y, θ) - lcf(X+Y, θ) )
    end
end

c(X::AbstractVector, Y::AbstractVector, θ::Real = 1) = (ecf(X+Y, θ) + ecf(X-Y, θ))/2 - ecf(X, θ)*ecf(Y, θ)

"""
    lcfAsymptDistr(X::AbstractVector, θ::Real = 1)

Asymptotic distribution of the empirical lcf.
"""
function lcfAsymptDistr(X::AbstractVector, θ::Real = 1)
    n = length(X)
    μ = lcf(X, θ)
    σ = 2/θ^2 * √c(X, X, θ) / ecf(X, θ)
    return Normal(μ, σ/√n)
end

"""
    codAsymptDistr(X::AbstractVector, Y::AbstractVector, θ::Real = 1, type::Symbol = :s)

Asymptotic distribution of the emprical codifference.
"""
function codAsymptDistr(X::AbstractVector, Y::AbstractVector, θ::Real = 1, type::Symbol = :s)
    type ∉ (:s, :+, :-) && error("Unrecognised type of codifference")
    n = length(X)
    μ = cod(X, Y, θ, type)
    σ = if type == :s
        ϕxmy, ϕxpy = ecf(X-Y, θ), ecf(X+Y, θ)
        1/2θ^2 * sqrt(c(X+Y, X+Y, θ)/ϕxpy^2 + c(X-Y, X-Y, θ)/ϕxmy^2 - 2c(X+Y, X-Y, θ)/(ϕxpy*ϕxmy))
    elseif type == :+
        ϕx, ϕy, ϕxmy = ecf(X, θ), ecf(X, θ), ecf(X-Y, θ)
        1/θ^2 * sqrt(c(X, X, θ)/ϕx^2 + c(Y, Y, θ)/ϕy^2 + c(X-Y, X-Y, θ)/ϕxmy^2 + 2c(X, Y, θ)/(ϕx*ϕy) - 2c(X, X-Y, θ)/(ϕx*ϕxmy) - 2c(Y, X-Y, θ)/(ϕy*ϕxmy))
    else
        ϕx, ϕy, ϕxpy = ecf(X, θ), ecf(X, θ), ecf(X+Y, θ)
        1/θ^2 * sqrt(c(X, X, θ)/ϕx^2 + c(Y, Y, θ)/ϕy^2 + c(X+Y, X+Y, θ)/ϕxpy^2 + 2c(X, Y, θ)/(ϕx*ϕy) - 2c(X, X+Y, θ)/(ϕx*ϕxpy) - 2c(Y, X+Y, θ)/(ϕy*ϕxpy))
    end
    return Normal(μ, σ/√n)
end


"""
    lcfConfInterval(X::AbstractVector, θ::Real = 1, p::Real = 0.95)

Asymptotic confidence interval of the empircal lcf.
"""
function lcfConfInterval(X::AbstractVector, θ::Real = 1, p::Real = 0.95)
    α = 1-p
    d = lcfAsymptDistr(X, θ)
    return (quantile(d, α/2), quantile(d, 1-α/2))
end

"""
    codConfInterval(X::AbstractVector, Y::AbstractVector, θ::Real = 1, p::Real = 0.95, type::Symbol = :s)

Asymptotic confidence interval of the empircal codifference.
"""
function codConfInterval(X::AbstractVector, Y::AbstractVector, θ::Real = 1, p::Real = 0.95, type::Symbol = :s)
    α = 1-p
    d = codAsymptDistr(X, Y, θ, type)
    return (quantile(d, α/2), quantile(d, 1-α/2))
end

end