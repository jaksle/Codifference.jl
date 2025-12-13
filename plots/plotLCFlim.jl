using Distributions, CairoMakie
using LaTeXStrings

##

p = 1/2
t = 1
fn1(w) = -2/t^2*log(p*exp(-t^2/2) + (1-p)*exp(-w*t))
fn2(w) = -2/t^2*log(p*exp(-t^2/2) + (1-p)*1/(1+w^2))
fn3(w) = -2/t^2*log(p*exp(-t^2/2) + (1-p)*exp(-t^2*w^2/2))
lim1 = -2/t^2*log(p*exp(-t^2/2) + (1-p)*1)
lim2 = 1 - 2*log(p)/t^2
ws = LinRange(0.01,100,10000)

with_theme(theme_latexfonts()) do
fig = Figure(size=(700,350))

ax = Axis(fig[1,1],
    xlabel = "x",
    ylabel = "pdf",
    title = "Bulk and tails of the pdf",
    limits = (-10,10,nothing,nothing),
)
xs = LinRange(-10,10,200)
p1 = 0.7pdf.(Normal(0,1),xs)
bulk = band!(xs,zeros(size(xs)), p1,
    alpha = 0.5,
    color = :blue,
)
p2 = 0.3pdf.(Normal(0,10),xs)
tail = band!(xs,zeros(size(xs)), p2,
    alpha = 0.5,
    color = :red,
)
total = lines!(xs,p1 .+ p2,
    color = :black,
    linestyle = :solid,
    linewidth = 0.7,
)
axislegend(ax,[PolyElement(color = :blue,alpha=0.25), PolyElement(color = :red,alpha=0.25), total],["bulk","tails", "mixture"],
    position = :rc,
)
ax = Axis(fig[1,2],
    xlabel = "w",
    ylabel = "lcf",
    xscale = log10,
    xminorgridvisible = true,
    xminorticksvisible = true,
    title = "The lcf as a function of the tail width",
    xminorticks = IntervalsBetween(5),
    limits = (0.01,100,nothing,nothing)
)
l1 = lines!(ax,ws,fn1.(ws),
    label = "Cauchy",
)
l2 = lines!(ax,ws,fn2.(ws),
    label = "Laplace",
)
l3 = lines!(ax,ws,fn3.(ws),
    label = "Gaussian"
)
h = hlines!(ax,[lim1,lim2],
    xmin = -2,
    linestyle = :dash,
    color = :black,
    label = "limit",
    alpha = 0.5
)

axislegend(ax,[l1,l2,l3],["Cauchy","Laplace","Gaussian"],
    position = :rc,
)

save("LCFlim.pdf",fig)
fig
end