using Distributions, CairoMakie


##



with_theme(theme_latexfonts()) do
    fig = Figure(size = (800, 400))
    ax = Axis(fig[1, 1], aspect = 1,
        xlabel = "x",
        ylabel = "y",
        #xgridvisible = false,
        #ygridvisible = false,
        xlabelsize= 20,
        ylabelsize = 20,
    )
    a, b = 0.02,1
    d1 = MvNormal([a 0; 0 b])
    d2 = MvNormal([b 0; 0 a])
    xs = LinRange(-2,2,200)
    z = [0.5*(pdf(d1,[x,y])+pdf(d2,[x,y])) for x in xs, y in xs]
    contour!(ax,xs,xs,z,
        levels = 10,
        color = :black,
        #linestyle = :dash,
        labels = false,

    )
    ax = Axis(fig[1, 2], aspect = 1,
        xlabel = "x",
        #xgridvisible = false,
        #ygridvisible = false,
        xlabelsize= 20,
        ylabelsize = 20,
    )

    c, s = 0.95,  1
    d1 = MvNormal([s s*c; s*c s])
    d2 = MvNormal([s -s*c; -s*c s])

    z = [0.5*(pdf(d1,[x,y])+pdf(d2,[x,y])) for x in xs, y in xs]
    contour!(ax,xs,xs,z,
        levels = 10,
        color = :black,
        #linestyle = :dash,
        labels = false,

    )
    save("cross.pdf",fig)
    fig
end

##
darkRed = colorant"#cc3434"

with_theme(theme_latexfonts()) do
    fig = Figure(size=(1200,400),figure_padding=(0,20,0,0))
    xlab = [L"10^{%$i}" for i in -3.4:0.2:-2.6]
    xlab[end÷2+1] = L"10^{-3}"
    xname = L"{D}\ [\mu m^2/s^\alpha]"
    ax = Axis(fig[1,1],
        xticks = (-3.4:0.2:-2.6,xlab),
        limits = (-3.4,-2.6, 0.4,1),
        title = "Step 1: points and predicted errors",
        titlesize = 20,
        xlabelsize= 20,
        ylabelsize = 20,
        xlabel = xname,
        ylabel = L"{\alpha}\ [1]",
        xgridvisible = false,
        ygridvisible = false,
    )
    gls = scatter!(ax,eB[1,:],eB[2,:],
        #markerstrokewidth=0,
        markersize=6,
        alpha = 0.15,
        color = :tomato, #darkRed,
        #label = "",

    )
    scatter!(ax, [log10(D0)],[2H0],
        marker='⨉',
        color=:black,
        markersize = 15,
    )

    C = sqrt(eM)
    θs = LinRange(0,2pi,200)
    xs =  @. C[1,1]*f(θs)+C[1,2]*g1(θs) + log10(D0)
    ys = @. C[2,1]*f(θs)+C[2,2]*g1(θs)+ 2H0
    conf = lines!(ax,xs,ys ,
        linewidth = 1.5,
        color = :black,
        linestyle = :dash,
        label = "",
    )


    ax2 = Axis(fig[1,2],
        xticks = (-3.4:0.2:-2.6,xlab),
        xlabel = xname,
        limits = (-3.4,-2.6, 0.4,1),
        title = "Step 2: density estimate",
        titlesize = 20,
        xlabelsize= 20,
        ylabelsize = 20,
    )
    heatmap!(ax2, den.x,den.y,den.density,
        colormap = :thermal,
    )
    scatter!(ax2, [log10(D0)],[2H0],
        marker='⨉',
        color=:black,
        markersize = 15,
    )

    ax3 = Axis(fig[1,3],
        xticks = (-3.4:0.2:-2.6,xlab),
        xlabel = xname,
        limits = (-3.4,-2.6, 0.4,1),
        title = "Step 3: deconvolution",
        titlesize = 20,
        xlabelsize= 20,
        ylabelsize = 20,
    )
    heatmap!(ax3, den.x,den.y,res)
    cross = scatter!(ax3, [log10(D0)],[2H0],
        marker='⨉',
        color=:black,
        markersize = 15,
    )

    axislegend(ax,[MarkerElement(color = :tomato, marker=:circle, alpha = 0.6, markersize = 12),  conf, cross,],["GLS","error 95% ellipse",L"exact ($D, \alpha$)"],
        position = :lt,
    )
    fig
    #save("decEx.pdf",fig)
end