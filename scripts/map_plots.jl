using PyCall
anim = pyimport("matplotlib.animation")
using DataFrames
import Interpolations: interpolate, extrapolate, BSpline, Constant, Flat
import ImageFiltering: imfilter, Kernel

import MapImages: nearest_pixel, nearest, grid, km_to_deg

### Gridding Functions for summing CSVs into bins ###
#
oob(row, col, out) = row < 1 || row > size(out, 1) || col < 1 || col > size(out, 2)

function bin_vals(demrsc, df; valcol=:sum15_17, loncol=:LongNAD27, latcol=:LatNAD27,
                  smoothing::Union{Nothing, <:Real}=nothing)
    out = zeros(size(demrsc))
    for (lon, lat, val) in eachrow(df[:, [loncol, latcol, valcol]])
        row, col = nearest_pixel(demrsc, lat, lon)
        oob(row, col, out) && continue
        out[row, col] += coalesce(val, 0)
    end
    
    return (!isnothing(smoothing) && smoothing > 0) ? imfilter(out, Kernel.gaussian(smoothing)) : out
end

function coarse_bin_vals(demrsc, df; step_km=nothing, digits=nothing, 
                         valcol=:sum15_17, loncol=:LongNAD27, latcol=:LatNAD27,
                         samesize=true, sum_vals=true)
    if !isnothing(step_km)
        lons, lats = coarse_grid_km(demrsc, step_km)
    elseif !isnothing(digits)
        lons, lats = coarse_grid_deg(demrsc, digits)
    else
        error("Need step_km or digits")
    end
    out = zeros(size(lats, 1), size(lons, 1))
    for (lon, lat, val) in eachrow(df[:, [loncol, latcol, valcol]])
        row = nearest(lats, lat)
        col = nearest(lons, lon)
        oob(row, col, out) && continue
        out[row, col] += (sum_vals ? val : 1)
    end
    return samesize ? resize_nearest(out, size(demrsc)...) : out
end

function resize_nearest(img, nrows, ncols)
    intp = interpolate(img, BSpline(Constant()))
    etpf = extrapolate(intp, Flat())
    return etpf(range(.5, stop=size(img, 1)+.5, length=nrows), range(.5, stop=size(img, 2)+.5, length=ncols))
end


"""Create a grid in DemRsc area with spacing `step_km`"""
function coarse_grid_km(demrsc, step_km::Union{AbstractFloat, Nothing}=(10/0.62))
    lons, lats = grid(demrsc, sparse=true)
    step = km_to_deg(step_km)
    return lons[1]:step:lons[end], lats[1]:(-step):lats[end]
end

"""Create a grid from a DemRsc on the 10^(-`digits`) grid lines"""
function coarse_grid_deg(demrsc, digits::Union{Int, Nothing}=2)
    step = 10.0^(-digits)
    lons, lats = grid(demrsc, sparse=true)
    lat1, lat2 = round.(extrema(lats), digits=digits)
    lon1, lon2 = round.(extrema(lons), digits=digits)
    return lon1:step:lon2, lat2:(-step):lat1
end

function plot_wells_per_mi(demrsc)
    oil_prod = CSV.read("oil_production.csv");
    oil_per_mi = coarse_bin_vals(demrsc, oil_prod, sum_vals=false, step_km=1/0.62);

    df_sub = df[!, coalesce.(df[:, :sum15_17], 0) .> 0]

    colors = broadcast(t -> t./256, [ (69, 117, 199, 256), (145, 191, 219, 256), (224, 243, 248, 256), (255, 255, 191, 255), (254, 224, 144, 256), (252, 141, 89, 256), (215, 48, 39, 256), ])
    plt.figure(); 
    plt.contourf(oil_per_mi, colors=colors, levels=[0, 1, 5, 10, 15, 25, 50, maximum(oil_per_mi)], origin="image", vmax=45, extend="max")
    plt.colorbar()
end

function sum_years(df, years...; outcol=:yearsum)
    df_out = copy(df)
    yearsyms = [Symbol(y) for y in years]
    df_out[!, outcol] = zeros(size(df_out, 1))

    for i = 1:length(yearsyms)
        df_out[!, outcol] += coalesce.(df[:, yearsyms[i]], 0)
    end
    return df_out
end


function plot_wells_per_mi(demrsc; years=[2015, 2016, 2017])
    oil_prod = CSV.read("oil_production.csv");

    df_sub = df[!, coalesce.(df[:, :sum15_17], 0) .> 0]

    oil_per_mi = coarse_bin_vals(demrsc, oil_prod, sum_vals=false, step_km=1/0.62);

    colors = broadcast(t -> t./256, [ (69, 117, 199, 256), (145, 191, 219, 256), (224, 243, 248, 256), (255, 255, 191, 255), (254, 224, 144, 256), (252, 141, 89, 256), (215, 48, 39, 256), ])
    plt.figure(); 
    plt.contourf(oil_per_mi, colors=colors, levels=[0, 1, 5, 10, 15, 25, 50, maximum(oil_per_mi)], origin="image", vmax=45, extend="max")
    plt.colorbar()
end

# E.g.
# m = rand(5, 5);
# stack = cat([m .+ .2i for i in 1:10]..., dims=3)
function animate_stack(stack; outname="test.gif", vm=maximum(stack), delay=200, cmap="seismic_wide_y",
                      geolist=nothing)

    fig, ax = plt.subplots()

    # function make_frame(i)
    #     plt.imshow(stack[:,:,i+1], vmax=vm, vmin=-vm)
    # end
    # myanim = anim.FuncAnimation(fig, make_frame, frames=size(stack,3), interval=20)
    # myanim[:save]("test2.mp4", bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])

    titles = !isnothing(geolist) ? string(geolist) : [string(i) for i in 1:size(stack,3)]
    for i = 1:size(stack, 3)
        ax.clear()
        ax.imshow(stack[:,:,i], vmax=vm, vmin=-vm, cmap=cmap)
        ax.set_title(titles[i])
        fname = @sprintf("testgif_%04d",i)
        println("Saving $(titles[i]) as $fname")
        fig.savefig(fname, bbox_inches="tight")
    end

    # Idk why, but you need to divide by 10 for it to really be ms for magick
    #run(`convert -delay $(interval/10) -loop 0 tmp0\*.png tmp.gif`)
    run(`magick -delay $(delay/10) -loop 0 testgif_\*.png $outname`)
    rm.(glob("testgif_*.png"))

    return outname
end

function animate_imgs_pts(df, date_list, stack::Union{MapImage, Nothing}=nothing;
                          loncol=:lon, latcol=:lat, sizecol=:mag, datecol=:datetime,
                          size_scale=4, outname="test.gif", alpha=.4, vm=6, c="r", 
                          delay=200, cmap="seismic_wide_y", extent=nothing, titles=nothing)
    # (-104.0, -103.001666673056, 30.901666673336, 31.90000000028)
    extent = isnothing(stack) ? extent : MapImages.grid_extent(stack)
    xmin, xmax, ymin, ymax = extent

    # fig, ax = plt.subplots()
    fig, ax = plot_states(extent)
    # titles = !isnothing(geolist) ? string(geolist) : [string(i) for i in 1:size(stack,3)]
    titles = isnothing(titles) ? string(date_list) : titles


    prevdate = Date(1,1,1)
    # for i=1:size(stack, 3)
    for i=1:length(date_list)
    # for i=1:4
        curdate = date_list[i]
        df_sub = df[(df[:, datecol] .< curdate) .& (df[:, datecol] .> prevdate), :]
        prevdate = curdate

        (lons, lats, sizes) = (df_sub[:, loncol], df_sub[:, latcol], df_sub[:, sizecol])
        sizes .*= size_scale

        ax.clear()
        fig, ax = plot_states(extent, fig, ax)

        ax.scatter(lons, lats, sizes, c=c, alpha=alpha, edgecolor="none", zorder=2.5, transform=ccrs.PlateCarree())
        !isnothing(stack) && ax.imshow(stack[:,:,i], vmax=vm, vmin=-vm, cmap=cmap, extent=extent)
        ax.set_title(titles[i])
        fname = @sprintf("testgif_%04d",i)
        println("Saving $(titles[i]) as $fname")
        fig.savefig(fname, bbox_inches="tight")
    end
    run(`magick -delay $(delay/10) -loop 0 testgif_\*.png $outname`)
    rm.(glob("testgif_*.png"))
    return fig, ax
end
# E.G.
# animate_imgs_pts(eqs_txar_big, Date(2001,1,1):Dates.Year(1):Date(2018,1,1), extent=(-105., -100., 29., 33.), delay=650, size_scale=8, titles=["TXAR earthquakes in $(string(year(d)))" for d in date_list])
#
# animate_imgs_pts(eqs_txar_big, Date(2001,1,1):Dates.Year(1):Date(2018,1,1), extent=MapImages.grid_extent(demrsc), delay=550, size_scale=8, titles=[string(year(d)) for d in date_list])
# animate_imgs_pts(eqs_txar, Date(2001,1,1):Dates.Year(1):Date(2017,1,1), extent=MapImages.grid_extent(demrsc), delay=350)
# animate_imgs_pts(eqs_txar, geolist, stack78)

# plt.figure(); plt.contourf(oil_per_mi, colors=colors, levels=[0, 1, 5, 10, 15, 25, 50, maximum(oil_per_mi)], origin="image", vmax=45, extend="max"); plt.colorbar()
ccrs = pyimport("cartopy.crs")
shpreader = pyimport("cartopy.io.shapereader")
sgeom = pyimport("shapely.geometry")
cfeature = pyimport("cartopy.feature")

function plot_states(extent, fig=nothing, ax=nothing; add_counties=true)
    fig = isnothing(fig) ? plt.figure() : fig
    ax = isnothing(ax) ? fig.add_axes([0, 0, 1, 1], projection=ccrs.LambertConformal()) : ax

    # extent = isnothing(stack) ? extent : MapImages.grid_extent(stack)
    ax.set_extent(extent, ccrs.Geodetic())

    shapename = "admin_1_states_provinces_lakes_shp"
    # states_shp = shpreader.natural_earth(resolution="110m", category="cultural", name=shapename)
    # states_shp = shpreader.natural_earth(resolution="50m", category="cultural", name=shapename)
    states_shp = shpreader.natural_earth(resolution="10m", category="cultural", name=shapename)

    # to get the effect of having just the states without a map "background"
    # turn off the outline and background patches
    ax.background_patch.set_visible(false)
    ax.outline_patch.set_visible(false)

    # records = shpreader.Reader(states_shp).records()
    # states = shpreader.Reader(states_shp).geometries()
    # ff = shpreader.Reader("/Users/scott/.local/share/cartopy/shapefiles/natural_earth/cultural/ne_10m_admin_1_states_provinces_lakes_shp.shp")
    ff = shpreader.Reader("/Users/scott/Downloads/cb_2018_us_state_5m/cb_2018_us_state_5m.shp")
    records = ff.records()
    states = ff.geometries()
    # records = ff.records()
    for (record, state) in zip(records, states)
    # for state in states
        # pick a default color for the land with a black outline,
        # this will change if the storm intersects with our track
        facecolor = record.attributes["NAME"] == "Texas" ? [0.9375, 0.9375, 0.859375] : [1, 1, 1, 0]
        # facecolor = [0.9375, 0.9375, 0.859375]
        edgecolor = "black"

        # if state.intersects(track)
        #     facecolor = "red"
        # elseif state.intersects(track_buffer):
        #     facecolor = "#FF7E00"
        # end

        ax.add_geometries([state], ccrs.PlateCarree(),
                          facecolor=facecolor, edgecolor=edgecolor)
    end
    pecos =(-103.49283 ,  31.4243)
    ax.plot(pecos..., "bo", ms=10, transform=ccrs.PlateCarree() ,zorder=3)
    plt.text((.09 .+ pecos)..., "Pecos", fontdict=Dict("size"=>14, "weight"=>"bold"), horizontalalignment="left", transform=ccrs.Geodetic())

    if add_counties
        reader = shpreader.Reader("./countyl010g_shp_nt00964/countyl010g.shp")
        counties = collect(reader.geometries())
        COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())
        ax.add_feature(COUNTIES, facecolor="none", edgecolor="gray")
    end

    return fig, ax
end

# oil_totals = CSV.read("PermianTotalByYear.csv"); oil_totals = oil_totals[1:end-1, :]
function animate_bar(xs, ys; outname="test.gif", delay=200, bg_color="b", active_color="r", titles=string.(xs),
                    bigfont=true)
    # years, values = oil_totals[:, :year], oil_totals[:, :oil_bbl]
    fig, ax = plt.subplots()
    if bigfont
        rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
        rcParams["font.size"] = 16
        rcParams["font.weight"] = "bold"
    end

    for i=1:length(xs)
        ax.clear()

        colors = fill(bg_color, length(xs)); colors[i] = active_color

        ax.bar(xs, values, color=colors)

        ax.set_xticks([2000, 2009, 2017])
        ax.set_yticks([0, 0.5, 1, 1.5, 2])
        ax.set_ylabel("Million bbl/day")


        ax.set_title(titles[i])
        fname = @sprintf("testgif_%04d",i)
        println("Saving $(titles[i]) as $fname")
        fig.savefig(fname, bbox_inches="tight")
    end
    run(`magick -delay $(delay/10) -loop 0 testgif_\*.png $outname`)
    rm.(glob("testgif_*.png"))
    return fig, ax
 end
