using Distributions, CairoMakie, LaTeXStrings

with_theme(theme_latexfonts()) do
    fig = Figure(size = (500, 420),
        #figure_padding = (10,10,10,10)
        )
    xs = LinRange(-2pi,2pi,200)
    ax = Axis(fig[1, 1],
        xlabel = "x",
        ylabel = "y",
        xticks = (-2pi:pi:2pi, [L"-2\pi", L"-\pi", L"0", L"\pi", L"2\pi"]),
        yticks = (-2pi:pi:2pi, [L"-2\pi", L"-\pi", L"0", L"\pi", L"2\pi"]),
        #xgridvisible = false,
        #ygridvisible = false,
        xlabelsize= 20,
        ylabelsize = 20,
    )
    hm = heatmap!(ax,xs,xs,sin.(xs) .* sin.(xs'))
    c, s = 0.90,  7
    d = MvNormal([s s*c; s*c s])
    z = [pdf(d,[x,y]) for x in xs, y in xs]
    contour!(ax,xs,xs,z,
        color=:white,
        linestyle = :dash,
        linewidth = 2,
        levels = 5,
    )
    Colorbar(fig[1, 2],limits = (-1, 1))
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    save("sinsin.pdf",fig)
    fig
end