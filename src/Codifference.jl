module Codifference

using StatsBase, Distributions

export ecf, lcf, cdf


"""
    ecf(X::AbstractVector, θ::Real)

Empirical characteristic function for random vectors with symmetric distribution. Returns real valued output
"""
ecf(X::AbstractVector, θ::Real) = mean( cos(θ*x) for x in X)

lcf(X::AbstractVector, θ::Real = 1) = θ ≈ 0 ? one(θ)*var(X) : -2/θ^2 * log(ecf(X,θ))

function cdf(X::AbstractVector, Y::AbstractVector, θ::Real = 1, type::Symbol = :s)
    if type == :s
        1/4 * ( lcf(X+Y, θ) - lcf(X-Y, θ) )
    elseif type == :+
        1/2 * ( lcf(X, θ) + lcf(Y, θ) - lcf(X-Y, θ) )
    elseif type == :-
        -1/2 * ( lcf(X, θ) + lcf(Y, θ) - lcf(X+Y, θ) )
    end
end

function cdfAsDistr(X::AbstractVector, Y::AbstractVector, θ::Real = 1, type::Symbol = :s)
    if type == :s
        1/4 * ( lcf(X+Y, θ) - lcf(X-Y, θ) )
    elseif type == :+
        1/2 * ( lcf(X, θ) + lcf(Y, θ) - lcf(X-Y, θ) )
    elseif type == :-
        -1/2 * ( lcf(X, θ) + lcf(Y, θ) - lcf(X+Y, θ) )
    end
end

end