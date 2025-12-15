using CairoMakie, Distributions, ProgressMeter
using .Codifference

## lcf

X = [rand(Exponential(1))*randn() for _ in 1:5*10^4]

S = zeros(10^4)

@showprogress for k in 1:10^4
    X = [rand(Exponential(1))*randn() for _ in 1:5*10^4]
    S[k] = lcf(X)
end

fig = Figure(); ax = Axis(fig[1,1])
d = lcfAsymptDistr(X)
d2 = Normal(mean(S),d.σ)
hist!(ax,S, normalization = :pdf)
xs = LinRange(minimum(S),maximum(S),200)
lines!(ax,xs, pdf.(d2,xs))
fig

## cdf
function gen(n,N)
    S = zeros(N)
    X0 = randn(2,10^5)
    X0 .= rand(Exponential(1),10^5)'
    @showprogress for k in 1:N
        X = randn(2,n)
        X .= rand(Exponential(1),n)'
        S[k] = cod(X[1,:],X[2,:],1, :s)
    end
    return X0, S
end

with_theme(theme_latexfonts()) do
fig = Figure(size = (800,400))

ax = Axis(fig[1,1],
    title = "Sample size = 100",
    xlabel = "codifference",
    ylabel = "pdf",
    limits = (0.3,2,nothing,nothing) 
)
n, N = 10^2, 10^4
X0, S = gen(n,N)
S2 = S[.!isnan.(S)]
d = codAsymptDistr(X0[1,:],X0[2,:],1, :s)
d2 = Normal(mean(S2),√1000*d.σ)
hist!(ax,S2, normalization = :pdf,
    bins = 50,
    label = "simulation"
)
xs = LinRange(d2.μ-4d2.σ,d2.μ+4d2.σ,200)
lines!(ax,xs, pdf.(d2,xs),
    label = "prediction",
    linestyle = :dash,
    color = :black,
)
axislegend(ax,
    #position = :lt,
    position = (0.9, 1.0)
)

ax = Axis(fig[1,2],
    xlabel = "codifference",
    title = "Sample size = 500" 
)
n, N = 500, 10^4
X0, S = gen(n,N)
S2 = S[.!isnan.(S)]
d = codAsymptDistr(X0[1,:],X0[2,:],1, :s)
d2 = Normal(mean(S2),√200*d.σ)
hist!(ax,S2, normalization = :pdf,
    bins = 50,
)
xs = LinRange(d2.μ-4d2.σ,d2.μ+4d2.σ,200)
lines!(ax,xs, pdf.(d2,xs),
    linestyle = :dash,
    color = :black,
)


ax = Axis(fig[1,3],
    xlabel = "codifference",
    title = "Sample size = 1000" 
)
n, N = 10^3, 10^4
X0, S = gen(n,N)
S2 = S[.!isnan.(S)]
d = codAsymptDistr(X0[1,:],X0[2,:],1, :s)
d2 = Normal(mean(S2),√100*d.σ)
hist!(ax,S2, normalization = :pdf,
    bins = 50,
)
xs = LinRange(d2.μ-4d2.σ,d2.μ+4d2.σ,200)
lines!(ax,xs, pdf.(d2,xs),
    linestyle = :dash,
    color = :black,
)
save("asymptDist.pdf",fig)
fig
end