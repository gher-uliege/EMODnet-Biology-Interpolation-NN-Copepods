module CopepodsNNplots 

using GeoArrays
using PolygonOps
using Downloads
using GeoDatasets
using Makie
using CairoMakie
using GeoMakie
using DIVAnd


"""
	create_geoaxis(fig)

Create a GeoAxis for the current figure `fig`
"""
function create_geoaxis(
    fig::Figure;
    proj = "+proj=merc",
    title = "",
    ii = 1,
    jj = 1,
    longridvalid = [],
    latgridvalid = [],
    domain = [-180.0, 180.0, -90.0, 90.0],
)
    ga = GeoAxis(
        fig[ii,jj],
        dest = proj,
        xgridstyle = :dashdot,
        xgridcolor = :grey,
        ygridstyle = :dashdot,
        ygridcolor = :grey,
        xticks = collect(-110:10:30),
        yticks = collect(10.0:5:80),
        xminorticks = longridvalid,
        yminorticks = latgridvalid,
        title = title,
        titlefont="Arial",
        limits = ((domain[1], domain[2]), (domain[3], domain[4])),
    )
    return ga::GeoAxis
end

"""
	decorate_map(ga, lonmask, latmask, mask)
"""
function decorate_map(ga)
    mask2 = copy(mask)
    mask2[mask2 .== 0] .= NaN;
    heatmap!(ga, lonmask, latmask, mask2, colormap = Reverse(:Greys_3), colorrange = (0, 3))
    contour!(ga, lonmask, latmask, mask, levels = [0], color = :black, linewidth = 0.25)
    xlims!(ga, domain[1], domain[2])
    ylims!(ga, domain[3], domain[4])
end

end #module