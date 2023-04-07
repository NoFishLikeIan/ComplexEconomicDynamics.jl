function plotvectorfield(g::Function, grid::Tuple; plotkwargs...)
    fig = plot()
    plotvectorfield!(fig, g, grid; plotkwargs...)
    return fig
end

"""
Plot vector field of a function g(x, y) = (u, v) on a grid of points (x, y) ∈ (xs, ys).
"""
function plotvectorfield!(figure, 
    g::Function, grid::Tuple; 
    rescale=1, plotkwargs...)

    xlims = extrema(xs)
    ylims = extrema(ys)

    N, M = length(xs), length(ys)
    xm = repeat(xs, outer=M)
    ym = repeat(ys, inner=N)

    field = g.(xm, ym)

    scale = rescale * (xlims[2] - xlims[1]) / min(N, M)
    u = @. scale * first(field)
    v = @. scale * last(field)

    steadystates = @. (u ≈ 0) * (v ≈ 0)

    u[steadystates] .= NaN
    v[steadystates] .= NaN

    z = norm.(field)

    quiver!(
        figure, xm, ym;
        quiver=(u, v), line_z=repeat(z, inner=4),
        aspect_ratio=1, xlims=xlims, ylims=ylims,
        c=:batlow, colorbar=false,
        plotkwargs...
    )

end

export plotvectorfield, plotvectorfield!