module Codifference

using StatsBase, Distributions

export ecf, lcf, cdf, lcfAsymptDistr, cdfAsymptDistr, lcfConfInterval, cdfConfInterval


"""
    ecf(X::AbstractVector, θ::Real)

Empirical characteristic function for random vectors with symmetric distribution. Returns real valued output
"""
ecf(X::AbstractVector, θ::Real) = mean( cos(θ*x) for x in X)

"""
    lcf(X::AbstractVector, θ::Real = 1)

Empirical log characteristic function of sample X.
"""
lcf(X::AbstractVector, θ::Real = 1) = θ ≈ 0 ? one(θ)*var(X) : -2/θ^2 * log(ecf(X,θ))

"""
    cdf(X::AbstractVector, Y::AbstractVector, θ::Real = 1, type::Symbol = :s)

Empirical codifference of sample X, Y.
"""
function cdf(X::AbstractVector, Y::AbstractVector, θ::Real = 1, type::Symbol = :s)
    if type == :s
        1/4 * ( lcf(X+Y, θ) - lcf(X-Y, θ) )
    elseif type == :+
        1/2 * ( lcf(X, θ) + lcf(Y, θ) - lcf(X-Y, θ) )
    elseif type == :-
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
    cdfAsymptDistr(X::AbstractVector, Y::AbstractVector, θ::Real = 1, type::Symbol = :s)

Asymptotic distribution of the emprical codifference.
"""
function cdfAsymptDistr(X::AbstractVector, Y::AbstractVector, θ::Real = 1, type::Symbol = :s)
    n = length(X)
    μ = cdf(X, Y, θ, type)
    σ = if type == :s
        ϕxmy, ϕxpy = ecf(X-Y, θ), ecf(X+Y, θ)
        1/2θ^2 * sqrt(c(X+Y, X+Y, θ)/ϕxpy^2 + c(X-Y, X-Y, θ)/ϕxmy^2 - 2c(X+Y, X-Y, θ)/(ϕxpy*ϕxmy))
    elseif type == :+
        ϕx, ϕy, ϕxmy = ecf(X, θ), ecf(X, θ), ecf(X-Y, θ)
        1/θ^2 * sqrt(c(X, X, θ)/ϕx^2 + c(Y, Y, θ)/ϕy^2 + c(X-Y, X-Y, θ)/ϕxmy^2 + 2c(X, Y, θ)/(ϕx*ϕy) - 2c(X, X-Y, θ)/(ϕx*ϕxmy) - 2c(Y, X-Y, θ)/(ϕy*ϕxmy))
    elseif type == :-
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
    cdfConfInterval(X::AbstractVector, Y::AbstractVector, θ::Real = 1, p::Real = 0.95, type::Symbol = :s)

Asymptotic confidence interval of the empircal codifference.
"""
function cdfConfInterval(X::AbstractVector, Y::AbstractVector, θ::Real = 1, p::Real = 0.95, type::Symbol = :s)
    α = 1-p
    d = cdfAsymptDistr(X, Y, θ, type)
    return (quantile(d, α/2), quantile(d, 1-α/2))
end

end