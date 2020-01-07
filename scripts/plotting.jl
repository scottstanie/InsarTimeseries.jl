using HDF5
import PyPlot
import Polynomials
import Dierckx: Spline2D
using PyCall
kml = pyimport("apertools.kml")
using Printf
using Glob

import MapImages
import MapImages: MapImage
import ImageFiltering: Kernel, imfilter

plt = PyPlot
include("./colors.jl")  # Custom cmaps
 
load_geolist_intlist = InsarTimeseries.load_geolist_intlist
baseline = InsarTimeseries.temporal_baseline

function save_paper_figure(fig, fname, axis_off=false, dpi=350)
    fig.tight_layout()
    axis_off && [a.set_axis_off() for a in fig.axes]
    println("Saving $fname")
    fig.savefig(fname, bbox_inches="tight", transparent=true, dpi=dpi)
end

function savefig_nomargin(fig, fname)
    NullLocator = pyimport("matplotlib.ticker")

    plt.subplots_adjust(0,0,1,1,0,0)
    for ax in fig.axes
        ax.axis("off")
        ax.margins(0,0)
        ax.xaxis.set_major_locator(NullLocator.NullLocator())
        ax.yaxis.set_major_locator(NullLocator.NullLocator())
    end
    fig.savefig(filepath, pad_inches=0, bbox_inches="tight")
end


function read_last(fname, dset, do_permute=true)
    sz = size(fname, dset)
    if length(sz) > 2 
        data = h5read(fname, dset, (:, :, sz[3]))[:, :, 1]
    else
        data = h5read(fname, dset)
    end
    return do_permute ? permutedims(data) : data
end

function plotsplit(fname; cmap="seismic_wide", vm=nothing, n=1, 
                   group="velos", title="", twoway=true, shift=0)
    vs = [read_last(fname, "$group/$ii") for ii in 1:n]
    # elseif group == "stack"
    #

    if isnothing(vm)
        vm = maximum([maximum(abs.(vv)) for vv in vs])
        println("Using $vm as max colorbar")
    end
    vmin, vmax = twoway ? (-vm, vm) : (0, vm)

    fig, axes = plt.subplots(1, n, squeeze=false, sharex=true, sharey=true)
    axim = nothing
    for ii = 1:n
        vi = vs[ii]
        axim = axes[ii].imshow(vi .+ shift, vmin=vmin, vmax=vmax, cmap=cmap)
    end
    fig.colorbar(axim, ax=axes[end])

    if isempty(title)
        title = "$fname: $group"
    end
    fig.suptitle(title)
    plt.show(block=false)
    return fig, axes, vs
end

function _get_vminmax(img, vm=nothing, vmin=nothing, vmax=nothing; twoway=true)
    img_nonan = filter(!isnan, img);
    vm = isnothing(vm) ? maximum(abs.(img_nonan)) : vm
    if twoway
        vmax = isnothing(vmax) ? vm : vmax
        vmin = isnothing(vmin) ? -vm : vmin
    else
        vmax = isnothing(vmax) ? vm : vmax
        vmin = isnothing(vmin) ? 0 : vmin
    end
    return vmin, vmax
end

function plot_img(m::MapImage, fig, ax; cmap="seismic_wide_y", vm=nothing, 
                  vmax=nothing, vmin=nothing, title="", shift=0,
                  twoway=true, use_lat=false, point_list::AbstractArray=[], label_list=[])
    vmin, vmax = _get_vminmax(m, vm, vmin, vmax, twoway=twoway)

    extent = use_lat ? MapImages.grid_extent(m) : nothing
    axim = ax.imshow(m .+ shift, vmin=vmin, vmax=vmax, cmap=cmap, extent=extent)
    fig.colorbar(axim, ax=ax)
    for (idx, point) in enumerate(point_list)
        label = length(label_list) >= idx ? label_list[idx] : nothing
        ax.plot(point..., "x", label=label, ms=10)
    end

    length(label_list) > 0 && fig.legend()

    ax.set_title(title)
    plt.show(block=false)
    return fig, ax, m
end

function plot_img(fname::AbstractString, dset::AbstractString="velos/1"; kwargs...) 
    fig, ax = plt.subplots()
    return plot_img(MapImage(fname, dset), fig, ax; title="$fname: $dset", kwargs...)
end

function plot_img_diff(f1, f2, d1="velos/1", d2="velos/1"; 
                       title1="$f1: $d1", title2="$f2: $d2",
                       vm1=nothing, vmd=4, kwargs...)
    m1, m2 =  MapImage(f1, d1), MapImage(f2, d2)
    return plot_img_diff(m1, m2; title1=title1, title2=title2, vm1=vm1, vmd=vmd, kwargs...), m1, m2
end
function plot_img_diff(a1::AbstractArray, a2::AbstractArray;
                       title1="", title2="",
                       vm1=nothing, vmd=4, kwargs...)
    fig, axes = plt.subplots(1, 3, sharex=true, sharey=true)
    plot_img(a1, fig, axes[1]; vm=vm1, title=title1, kwargs...)
    plot_img(a2, fig, axes[2]; vm=vm1, title=title2, kwargs...)
    # Always do twoway for diff
    vmin, vmax = _get_vminmax(a1 - a2, vmd, twoway=true)
    plot_img(a1 - a2, fig, axes[3]; title="diff (1 - 2)", vmin=vmin, vmax=vmax, kwargs...)
    return fig, axes
end

function plot_regs(pixel, geolist, intlist; alpha=100, title="", L1=false)
    plt.figure()
    colors = ["b", "g", "c", "r"]
    prunes = [true, true, false, false]
    alphas = [alpha, 0.0, alpha, 0.0]
    labels = ["prune/$alpha", "prune/0", "no prune/$alpha", "no prune/0.0"]
    constant = false
    for ii = 1:4
        soln_velos, igram_count, geo_clean = InsarTimeseries.calc_soln(pixel, geolist, intlist, 1.0, alphas[ii],
                                                                       constant, L1=L1, prune=prunes[ii])
        phi_arr = InsarTimeseries._unreg_to_cm(soln_velos, geo_clean, geolist);
        plt.plot(geolist, phi_arr, colors[ii]*"-x", label=labels[ii])
    end
    plt.title(title)
    plt.legend()
    plt.show(block=false)
end

function _get_pixel_data(rowcol::Tuple{Int, Int}, unwfile="unw_stack.h5", max_temp=800, max_date=nothing)
    geolist, intlist, igram_idxs = load_geolist_intlist(unwfile, "geolist_ignore.txt", max_temp, max_date=max_date)
    B = InsarTimeseries.build_B_matrix(geolist, intlist);
    Blin = sum(B, dims=2);
    unw_vals = get_stack_vals(unwfile, rowcol..., 1, "stack_flat_shifted", igram_idxs)
    return geolist, intlist, igram_idxs, Blin, unw_vals
end

function plot_grouped_by_day(geo, int, vals, nsigma=0, min_spread=2)
    means = InsarTimeseries.mean_abs_val(geo, int, vals)
    println("median val: $(median(means))")
    fig, ax = plt.subplots()
    ax.scatter(geo, means, label="means by date")
    if nsigma > 0
        low, high = InsarTimeseries.two_way_cutoff(means, nsigma)
        ax.plot(geo, ones(length(geo)) * median(means), label="median")
        ax.plot(geo, ones(length(geo)) * high, label="$nsigma sigma MAD cutoff")
        ax.plot(geo, ones(length(geo)) * (median(means) + min_spread), label="min. spread")
        ax.plot(geo, ones(length(geo)) * nsigma*std(means), label="$nsigma sigma stddev cutoff")
    end
    fig.legend()
end

function plot_grouped_by_day(rowcol::Tuple{Int, Int}, unwfile="unw_stack.h5"; 
                             max_temp=800, max_date=nothing, nsigma=3, min_spread=0)
    geolist, intlist, igram_idxs, _, unw_vals = _get_pixel_data(rowcol, unwfile, max_temp, max_date)
    plot_grouped_by_day(geolist, intlist, unw_vals, nsigma, min_spread)
end


function plot_big_days(geolist, intlist, vals; to_cm=true, label=nothing, nsigma=4, 
                       min_spread=0, color=nothing, legend=true)

    B = sum(InsarTimeseries.prepB(geolist, intlist), dims=2);
    # means = InsarTimeseries.mean_abs_val(geolist, intlist, v);
    # high_days = sort(collect(zip(means, geolist)), rev=true)[1:5]
    high_days = InsarTimeseries.nsigma_days(geolist, intlist, vals, nsigma, min_spread)
    return plot_days(geolist, intlist, vals, B, high_days, to_cm=to_cm, label=label, color=color, legend=legend)
end

function plot_big_days(rowcol::Tuple{Int, Int}, unwfile="unw_stack.h5"; max_temp=800, max_date=nothing,
                       to_cm=true, label=nothing, nsigma=4, min_spread=0, color=nothing, legend=true)
    
    geolist, intlist, igram_idxs, _, unw_vals = _get_pixel_data(rowcol, unwfile, max_temp, max_date)
    plot_big_days(geolist, intlist, unw_vals; to_cm=to_cm, 
                  label=label, nsigma=nsigma, min_spread=min_spread, color=color)
end

function plot_days(geolist, intlist, vals, B, date_arr; to_cm=true, label=nothing, color=nothing, legend=true)
    scale = to_cm ? InsarTimeseries.PHASE_TO_CM : 1
    v = vals .* scale   

    fig, ax = plt.subplots()
    ax.scatter(B, v, label=label)

    ym = maximum(abs.(v)) * 1.2
    ylims = to_cm ? (-ym, ym) : (0, ym)  # correlation: 0 as bottom

    # for ii in 1:3
    for (ii, dd) in enumerate(date_arr)
        # hh = highest_days[ii][2]
        ax.scatter(InsarTimeseries.Blins_by_date(dd, intlist, B),
                   InsarTimeseries.vals_by_date(dd, intlist, v),
                   color=color,
                   label="date=$dd")
        ax.set_ylim(ylims)
    end
    legend && fig.legend()
    plt.show(block=false)
    return fig, ax
end


# Plot gps east up data
function plot_eu(station_name, insar_slopes=nothing; marker=".", title=station_name, ylim=(-2.5, 2.5))
    # Defined in find_abs_shift TODO fix this
    dts, east, north, up = get_gps_enu(station_name)

    insar_east, insar_up = !isnothing(insar_slopes) ? insar_slopes : (nothing, nothing)

    xticks = [dts[1], dts[Int(round(end//2))], dts[end]]

    fig, axes = plt.subplots(1, 2)
    axes[1].plot(dts, east, marker, label="gps")
    axes[1].set_title("East")
    axes[1].set_ylabel("cm")
    axes[1].set_ylim(ylim)
    axes[1].set_xticks(xticks)
    if !isnothing(insar_east)
        _plot_insar_line(axes[1], dts, insar_east)
    end

    axes[2].plot(dts, up, marker, label="gps")
    axes[2].set_title("Up")
    axes[2].set_ylim(ylim)
    axes[2].set_xticks(xticks)
    if !isnothing(insar_up)
        _plot_insar_line(axes[2], dts, insar_up)
    end

    fig.suptitle(title)
    return fig, axes
end
plot_gps_eu = plot_eu

function plot_enu(station_name, marker=".", title=station_name, ylim=(-2.5, 2.5))
    # Defined in find_abs_shift TODO fix this
    dts, east, north, up = get_gps_enu(station_name)
    xticks = [dts[1], dts[Int(round(end//2))], dts[end]]

    fig, axes = plt.subplots(3, 1)
    labels = ["East", "North", "Up"]

    for (ax, label, data) in zip(axes, labels, (east, north, up))
        ax.plot(dts, data, marker)
        ax.set_ylabel("$label (cm)")
        ax.set_ylim(ylim)
        ax.set_xticks(xticks)
        if ax != axes[end]
            plt.setp(ax.get_xticklabels(), visible=false)
        end
    end

    fig.suptitle(title)
    return fig, axes
end
plot_gps_enu = plot_enu

function _plot_insar_line(ax, dts, insar_slope; label="insar")
    # Offset to be zero centered for 3 years, slope to be per day
    p = Polynomials.Poly([-3(insar_slope / 2), insar_slope / 365])
    dd = InsarTimeseries._get_day_nums(dts)
    vals = p(dd)
    ax.plot(dts, vals, label=label)
end

function plot_all(east, up)
    station_overlap = ["TXMH", "TXFS", "TXAD", "TXS3", "NMHB"]
    slopes = eastups(station_overlap, east, up) ./ 10
    slope_by_station = collect(zip(slopes...))
    for (idx, name) in enumerate(station_overlap)
        plot_eu(name, slope_by_station[idx])
    end
end


# TODO: move these somehere
function count_dates(fname, dset)
    geolist = Sario.load_geolist_from_h5(fname, dset)
    excls = Sario.load(fname, dset_name="excluded/1")
    return geolist, excls
end 

function get_excl_hist(excl_map, geolist; as_pct=true, nonmasked_pixels=length(excl_map))
    counts = zeros(length(geolist))
    for idx in eachindex(excl_map)
        was_excluded = .!isnothing.(indexin(geolist, excl_map[idx]))
        counts .+= was_excluded
    end
    return as_pct ? counts ./ nonmasked_pixels : counts
end

function hist_change_from_outliers(station_name_list, fname, dset="velos/1";
                                   vm=20, bins::Int=70, ref_station="TXKM", maxdays=800,
                                  nrows=2)
    # geolist, excls = count_dates(fname, dset)
    # @time excl_map = InsarTimeseries.decode_excluded(excls, geolist);
    p2c = InsarTimeseries.PHASE_TO_CM

    plot_stations = [s for s in station_name_list if s != ref_station]
    nplots = length(plot_stations)  # Not plotting reference station

    # plt.xlabel("Centimeters Predicted by InSAR")
    # plt.ylabel("Interferogram count")
    # plt.title("InSAR Noise at Station NMHB")

    # m, n = 3, 5
    m, n = nrows, Int(ceil(nplots/nrows))
    # fig, axes = plt.subplots(m, n)
    fig, axes = plt.subplots(1, 1, squeeze=false)
    # for ii in 1:m
        # for jj in 1:n
            # idx = (ii-1)*n + jj
            # idx > length(station_name_list) && continue
    for (ax, name) in zip(axes, plot_stations)

        #$ name = plot_stations[idx]
        name == ref_station && continue

        # geolist, intlist, valid_igram_indices = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", maxdays)
        # unw_vals = get_stack_vals("unw_stack.h5", name, 1, "stack_flat_shifted", valid_igram_indices)
        unw_vals = h5read("gps_pixels_78.h5", name)
        geolist = Sario.load_geolist_from_h5("gps_pixels_78.h5")
        intlist = Sario.load_intlist_from_h5("gps_pixels_78.h5")

        g2, i2, u2 = InsarTimeseries.remove_outliers(geolist, intlist, unw_vals)

        # ax = axes[ii, jj]
        # ax = axes[idx]

        n1, _, _ = ax.hist(unw_vals .* p2c, range=(-vm, vm), bins=bins, label="All data")
        n1, _, _ = ax.hist(u2 .* p2c, range=(-vm, vm), bins=bins, label="Outliers removed")
        ax.set_xlabel("Centimeters")
        # ax.set_title("$name")
        ax.set_title("Interferogram count for $name")
        # ax.grid("on")
        println("$name: Before std. dev : $(std(unw_vals)), after: $(std(u2))")
    end
    # plt.legend()
    plt.show(block=false)
end

function plot_solution_grouped(unw_stack; vm=30, cmap="seismic_wide", baselines=200:100:600)
    solns = image_solution_grouped(unw_stack)
    n = length(solns) 
    fig, axes = plt.subplots(1, n)
    axim = nothing
    for i = 1:n
        axim = axes[i].imshow(solns[i], vmin=-vm, vmax=vm, cmap=cmap)
        axes[i].set_title("Max days = $(baselines[i])")
    end
    # plt.matplotlib.colorbar.make_axes([ax for ax in axes.flat])
    fig.colorbar(axim, ax=axes, shrink=0.75)
    return fig, axes
end


function image_solution_grouped(unw_stack; baselines=200:100:600)
    geolist, intlist, valid_igram_indices= load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 2000)
    unw_valid = unw_stack[:, :, valid_igram_indices]
    # Ball = InsarTimeseries.build_B_matrix(geolist, intlist);
    # Blinall = sum(Ball, dims=2);
 
    solns = []
    for b2 in baselines
    # for (b1, b2) in zip(baselines[1:end-1], baselines[2:end])
        @show b2
        # idxs = (baseline(intlist) .> b1) .& (baseline(intlist) .< b2)
        idxs = (baseline(intlist) .< b2)
        cur_ints = intlist[idxs]
        cur_B = sum(InsarTimeseries.build_B_matrix(geolist, cur_ints), dims=2)
        layers = @view unw_valid[:, :, idxs]
        # push!(solns, @view median(layers, dims=3)[:, :, 1])
        push!(solns, _calcstack(layers, intlist[idxs]))
    end
    return solns
end

function _calcstack(unw_stack, igrams)
    timediffs = baseline(igrams)
    phase_sum = sum(unw_stack, dims=3)[:, :, 1]
    avg_velo = phase_sum ./ sum(timediffs)
    # Finally, save as a mm/year velocity
    return InsarTimeseries.P2MM * avg_velo
end

function image_line(img, rowcol1, rowcol2)
    n = Int(round(sqrt(sum((rowcol2 .- rowcol1).^2 ))))
    # Note on Dierckx: "it is required that size(z) == (length(x), length(y))"
    XX = 1:size(img, 1)
    YY = 1:size(img, 2)
    intp = Spline2D(XX, YY, img)

    r1, c1 = rowcol1
    r2, c2 = rowcol2
    xs = range(r1, stop=r2, length=n)
    ys = range(c1, stop=c2, length=n)
    return intp(xs, ys)
end

function plot_image_line(img, rowcol1, rowcol2; plotkwargs...)
    line = image_line(img, rowcol1, rowcol2)
    fig, axes = plt.subplots(1, 2)
    axes[1].imshow(img; plotkwargs...)
    r1, c1 = rowcol1
    r2, c2 = rowcol2
    axes[1].plot(c1, r1, "gx", ms=10)
    axes[1].plot(c2, r2, "rx", ms=10)
    axes[1].plot([c1, c2], [r1, r2], "k")

    axes[2].plot(line)
end

# TODO: extract if not MapImage
function plot_image_line(img::MapImage, rowcol1, rowcol2; plotkwargs...)
    line = image_line(img, rowcol1, rowcol2)
    km_line_dist = MapImages.rowcol_to_dist(img, rowcol1, rowcol2)

    fig, axes = plt.subplots(1, 2)
    axes[1].imshow(img; plotkwargs...)
    r1, c1 = rowcol1
    r2, c2 = rowcol2

    axes[1].plot(c1, r1, "gx", ms=10)
    axes[1].plot(c2, r2, "rx", ms=10)
    axes[1].plot([c1, c2], [r1, r2], "k")

    axes[2].plot(0, line[1], "gx", ms=10)
    axes[2].plot(km_line_dist, line[end], "rx", ms=10)
    axes[2].plot(range(0, stop=km_line_dist, length=length(line)), line, "b.--")
    axes[2].set_ylabel("cm")
    axes[2].set_xlabel("km")
end

function plot_gps_los(name, insar_mm=nothing; ref="TXKM", start_date=Date(2014,11,1), 
                      end_date=(2019,2,1), ylim=(-3, 3), title="", bigfont=false, offset=true,
                      lw=5, gps_color="#86b251")
    fig, ax = plt.subplots()
    dts, los = get_gps_los(name, reference_station=ref,
                           start_date=start_date, end_date=end_date)
    day_nums = _get_day_nums(dts)

    if bigfont
        rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
        rcParams["font.size"] = 20
        rcParams["font.weight"] = "bold"
    end

    ax.plot(dts, los, ".", color=gps_color, markersize=7, label="GPS")

    if !isnothing(insar_mm)
        insar_cm_day = insar_mm / 365 / 10
        full_defo = insar_cm_day * (dts[end] - dts[1]).value
        bias = offset ? -full_defo/2 : 0
        ax.plot(dts, bias .+ day_nums .* insar_cm_day, "r", lw=lw, label="Insar")
    end

    # Or if you want different settings for the grids:
    ax.grid(which="major", alpha=0.5)
    ax.set_xticks(dts[1]:Dates.Day(365):dts[end])
    ax.set_yticks([-2, 0, 2])

    ax.set_ylim(ylim)
#     for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#              ax.get_xticklabels() + ax.get_yticklabels()):
#     item.set_fontsize(20)
    return fig, ax
end

function plot_los_maps(los_map_file="los_map.h5"; cmap="jet")
    stack = MapImages.MapImage(los_map_file, "stack")
    fig, axes = plt.subplots(1, 3, sharex=true, sharey=true)

    titles=["east", "north", "up"]
    for i = 1:3
        s = stack[:,:,i]
        axim = axes[i].imshow(abs.(s), cmap=cmap, vmin=0, vmax=1)
        fig.colorbar(axim, ax=axes[i])
        t = titles[i]
        axes[i].set_title("$t: extrema: $(round.(extrema(s), digits=2))")
    end
    return fig, axes, stack
end

function plot_15_vs_18(fnames=["velocities_2016_linear_max700.h5", "velocities_current.h5", "velocities_2018_linear_max700.h5"];
                       dset="velos_shifted/1", vm=6, cmap="seismic_wide_y", outnames=[])
    @time include("/home/scott/repos/MapImages/src/plotting.jl")
    fig, axes = plt.subplots(1, 3, sharex=true, sharey=true)
    for (idx, f) in enumerate(fnames)
        m = MapImage(f, dset)
        g = Sario.load_geolist_from_h5(f, dset)
        days = (g[end] - g[1]).value
        cum_map = m ./ 3650 * days
        axes[idx].imshow(cum_map, cmap=cmap, vmin=-vm, vmax=vm, extent=MapImages.grid_extent(m))
        if !isempty(outnames)
            fo = outnames[idx]
            plt.imsave(fo, cum_map, cmap=cmap, vmin=-vm, vmax=vm, format=strip(Sario.get_file_ext(fo), '.'), dpi=400)
        end
    end
    return fig, axes
end

function save_img_figure(outfile, array, levels, colors)
    ext = Sario.get_file_ext(outfile)
    @show outfile
    cmap_, norm_ = plt.matplotlib.colors.from_levels_and_colors(levels=levels, colors=colors)
    plt.imsave(outfile, norm_(array), cmap=cmap_, format=strip(ext, '.'))
    # plt.imsave(outfile, array, cmap=cmap, vmin=vmin, vmax=vmax, format=ext.strip('.'))
end


function save_img_geotiff(outfile::AbstractString, img, demrsc; cmap="seismic_wide_y", vm=nothing)
    ext = Sario.get_file_ext(outfile)
    vmin, vmax = _get_vminmax(img, vm)
    @show outfile
    outim = Float64.(img)
    outim[outim .== 0] .= NaN
    plt.imsave("tmp.png", outim, cmap=cmap, vmin=vmin, vmax=vmax, dpi=400, format="png")
    kml.create_geotiff(rsc_data=Dict(demrsc), img_filename="tmp.png", outfile=outfile)
    # rm("tmp.png")
end

save_img_geotiff(f, m::MapImage; kwargs...) = save_img_geotiff(f, m.image, m.demrsc; kwargs...)
save_img_geotiff(f, m::MapImage, lats::Tuple{AbstractFloat, AbstractFloat},
                 lons::Tuple{AbstractFloat, AbstractFloat}; kwargs...) = save_img_geotiff(f, m[lats, lons].image, m[lats, lons].demrsc; kwargs...)


function plot_stack_ts(stack, geolist, row, col; cmap="seismic_wide_y", vm=6, plotkwargs...)
    fig, axes = plt.subplots(1, 2)
    axim = axes[1].imshow(stack[:, :, end]; cmap=cmap, vmax=vm, vmin=-vm, plotkwargs...)
    fig.colorbar(axim, ax=axes[1])
    # axes[1].plot(col, row, "kx", ms=10)
    axes[1].plot(row, col, "kx", ms=10)
    plot_stack_ts(fig, axes[2], stack, geolist, row, col)
    return fig, ax
end

function plot_stack_ts(fig, ax, stack, geolist, row, col; label=nothing, title="")
    label = isnothing(label) ? "$row, $col" : label
    ax.plot(geolist, stack[row, col, :], "kx-", label=label)
    ax.set_ylabel("cm")
    ax.set_title(title)
    return fig, ax
end

extremanan(x) = extrema(replace(x, NaN => 0))

_num_days(g) = (g[end] - g[1]).value

import StatsBase: countmap
import Images: label_components
_max_label(labels) = sort([tup for tup in countmap(vec(labels)) if tup[1] != 0], by=x->x[2], rev=true)[1][1]

function save_pnas_images(;years=[2016, 2017, 2018], 
                          fnames=["velocities_$(year)_current.h5" for year in years],
                          # fnames=["velocities_$(year)_linear_max800_sigma4.h5" for year in years],
                          outdir=".",
                          # zoomlats=(31.6, 30.9), zoomlons=(-103.8, -103.),  # Zoomed just to stripes
                          zoomlats=(31.9, 31.0), zoomlons=(-103.8, -102.85),  # Covers txkm/txmh
                          latbounds=(32.8, 30.7041),
                          dset="velos_shifted/1", cmap1="seismic_wide_y", cmap2="seismic_wide_y", # cmap2=rdylbl, 
                          vm1=10, vm2=10,
                          dem_mask_height=nothing)
    maxdays = maximum([_num_days(Sario.load_geolist_from_h5(f, dset)) for f in fnames])
    @show maxdays
    vm1 *= (maxdays / 3650)
    vm2 *= (maxdays / 3650)
    @show vm1, vm2 

    m1, m2 = nothing, nothing
    for f in fnames
        @show f
        days = _num_days(Sario.load_geolist_from_h5(f, dset))
        @show days

        m1 = MapImage(f, dset)
        m1[m1 .== 0] .= NaN    

        if !isnothing(dem_mask_height)
            dem = Sario.load("elevation_looked.dem")
            labels = label_components(dem .< 1200)
            maxl = _max_label(labels)
            m1[labels .!= maxl] .= NaN
        end

        m1 = m1[latbounds, (-110., -90.)] # no lon bounds really
        @show MapImages.grid_extent(m1)

        @show extremanan(m1)
        m1 .*= (days / 3650)  # Cumulative, in cm
        @show extremanan(m1)
        m2 = m1[zoomlats, zoomlons]
        out1 = joinpath(outdir, replace(f, ".h5"  => ".tif"))
        # out2 = replace(f, ".h5"  => "_inset.tif")
        out2 = joinpath(outdir, replace(f, ".h5"  => "_inset_gps.tif"))

        # Save as both png and tif
        save_img_geotiff(out1, m1, cmap=cmap1, vm=vm1)
        # save_img_geotiff(out2, m2, cmap=cmap2, vm=vm2)
        plt.imsave(replace(out1, ".tif" => ".png"), m1, cmap=cmap2, vmin=-vm2, vmax=vm2, dpi=150)
        plt.imsave(replace(out2, ".tif" => ".png"), m2, cmap=cmap2, vmin=-vm2, vmax=vm2, dpi=150)
    end
    # Now save one with a colorbar
    fig, axes = plt.subplots(1, 2)
    axim = axes[1].imshow(m1, cmap=cmap1, vmax=vm1, vmin=-vm1)
    fig.colorbar(axim, ax=axes[1])
    axim = axes[2].imshow(m2, cmap=cmap2, vmax=vm2, vmin=-vm2)
    fig.colorbar(axim, ax=axes[2])

    cbar_name = "colorbars.png"
    save_paper_figure(fig, cbar_name, false, 400)
end


function plot_3d(img::MapImage; smooth=true, vm=nothing, colorbar=true)
    plt.using3D()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    XX, YY = MapImages.grid(img)
    img = smooth ?  imfilter(img, Kernel.gaussian(4)) : img
    vm = isnothing(vm) ? maximum(abs.(filter(!isnan, img))) : vm

    surf = ax.plot_surface(XX, YY, img.image, cmap="coolwarm", linewidth=0, 
                           antialiased=true, vmin=-vm, vmax=vm)

    ax.set_zlim(-vm, vm)

    colorbar && fig.colorbar(surf, shrink=0.5, aspect=5)
    return fig, ax
end

function pnas_outlier_figure(station_name="TXMC", outlier_color="r", h=3, w=3.5, bins=40)
    # Formatting functions:
    _remove_ticks(ax) = ax.tick_params(axis="both", which="both", bottom=false, top=false, labelbottom=false,
                                       left=false, right=false, labelleft=false)
    _set_figsize(fig, h=h, w=w) = (fig.set_figheight(h), fig.set_figwidth(w))

    # unw_vals = [h5read("gps_pixels.h5", name) for name in station_name_list78];
    unw_vals = h5read("gps_pixels_78.h5", station_name)
    geolist = Sario.load_geolist_from_h5("gps_pixels_78.h5");
    intlist = Sario.load_intlist_from_h5("gps_pixels_78.h5");

    end_date = Date(2018,1,1)
    geolist18 = geolist[geolist .> end_date]
    idxs = InsarTimeseries._good_idxs(geolist18, intlist)
    intlist17 = intlist[idxs]

    geolist17 = geolist[geolist .< end_date]
    unw_vals17 = unw_vals[idxs]
    
    # 1. plot grouped by date
    fig, ax = plot_big_days(geolist17, intlist17, unw_vals17, legend=false, color=outlier_color)
    _remove_ticks(ax)
    _set_figsize(fig)
    fig.savefig("outliers2.pdf", bbox_inches="tight", transparent=true, dpi=100)

    # 2. cm vs baseline igram outlier plot
    # function plot_grouped_by_day(geo, int, vals, nsigma=0, min_spread=2)
    geo = geolist17
    means = InsarTimeseries.mean_abs_val(geo, intlist17, p2c * unw_vals17)
    println("median val: $(median(means))")
    fig, ax = plt.subplots()
    ax.scatter(geo, means, label="means by date")
    nsigma = 4
    low, high = InsarTimeseries.two_way_cutoff(means, nsigma)
    ax.plot(geo, ones(length(geo)) * high, label="$nsigma sigma MAD cutoff")
    bigidx = means .> high 
    ax.scatter(geo[bigidx], means[bigidx], c=outlier_color)
    _remove_ticks(ax)
    _set_figsize(fig)  # Double the width here
    fig.savefig("outliers1.pdf", bbox_inches="tight", transparent=true, dpi=100)


    # 3. line  plot of insar after outliers vs GPS
    Blin = sum(InsarTimeseries.prepB(geolist17, intlist17), dims=2);
    before = (Blin \ unw_vals17) * p2c * 365 * 10  # Change to mm/year

    start_date=Date(2015,1,1)
    ref = nothing
    lw = 5
    fig, ax = plot_gps_los(station_name, before; ref=ref, lw=lw,
                           end_date=end_date, start_date=start_date,
                           ylim=(-3, 3), title="", bigfont=false, offset=false)

    # Now solve without outliers
    dts, los = get_gps_los(station_name, reference_station=ref, end_date=end_date, start_date=start_date)
    day_nums = _get_day_nums(dts)
    idxs = InsarTimeseries._good_idxs(geolist17[bigidx], intlist17)
    vals_after_outliers = unw_vals17[idxs]
    after = p2c * (Blin[idxs] \ vals_after_outliers)  # cm/day
    full_defo = after * (dts[end] - dts[1]).value
    # ax.plot(dts, -full_defo/2 .+ day_nums .* after, "c", lw=3, label="Insar")
    ax.plot(dts, day_nums .* after, "C0", lw=lw)
    _remove_ticks(ax)
    _set_figsize(fig, h, 2.3w)
    fig.savefig("outliers3.pdf", bbox_inches="tight", transparent=true, dpi=100)

    println("Before: $before mm/yr, $(before * (dts[end] - dts[1]).value /3650) total defo")
    println("After: $(after * 3650) mm/yr, $full_defo total defo")

    # 4. histogram before/after outlier
    fig, ax = plt.subplots()

    # ax = axes[ii, jj]
    # ax = axes[idx]
    vm = maximum(abs.(unw_vals17 * p2c))
    n1, _, _ = ax.hist(unw_vals17.* p2c, range=(-vm, vm), bins=bins, label="All data", color="red")
    n1, _, _ = ax.hist(vals_after_outliers .* p2c, range=(-vm, vm), bins=bins, label="Outliers removed", color="C0")
    _remove_ticks(ax)
    _set_figsize(fig)
    fig.savefig("outliers4.pdf", bbox_inches="tight", transparent=true, dpi=100)




    stds78 = Float64.(std.([h5read("gps_pixels_78.h5", name) for name in station_name_list78]))
    dists78 = [MapImages.latlon_to_dist(MapImages.station_latlon("TXKM"), MapImages.station_latlon(n)) for n in station_name_list78]
    df78 = DataFrame(dists=dists78, stds=stds78, names=station_name_list78)

    stds85 = Float64.(std.([h5read("gps_pixels_85.h5", name) for name in station_name_list85]))
    dists85 = [MapImages.latlon_to_dist(MapImages.station_latlon("TXKM"), MapImages.station_latlon(n)) for n in station_name_list85]
    df85 = DataFrame(dists=dists85, stds=stds85, names=station_name_list85)

    sm78 = lm(@formula(stds ~ dists), df78)
    a1, b1 = coef(sm78)
    sm85 = lm(@formula(stds ~ dists), df85)
    a2, b2 = coef(sm85)

    fig, ax = plt.subplots()
    ax.scatter(dists78, stds78, s=40)
    xs = 1:maximum(dists78)
    ax.plot(xs, a1 .+ b1 .* xs, lw=3)

    ax.scatter(dists85, stds85, s=40)
    xs = 1:maximum(dists85)
    ax.plot(xs, a2 .+ b2 .* xs, lw=3)

    ax.set_xlabel("km from TXKM")
    ax.set_ylabel("Std. Dev. of unw values")
    _remove_ticks(ax)
    _set_figsize(fig)
    fig.savefig("outliers5.pdf", bbox_inches="tight", transparent=true, dpi=100)

    return (geolist17, intlist17, unw_vals17)
end

_get_day_nums(dts::AbstractArray{Date, 1}) = [( d - dts[1]).value for d in dts]


function pnas_gps_error_table()
    p78 = "/data1/scott/pecos/path78-bbox2/igrams_looked/"
    p85 = "/data4/scott/path85/stitched/igrams_looked/"

   
    errors = []
    # errors85_4 = get_gps_error(p85*"velocities_2018_current.h5", all_stations, dset="velos_shifted/1", verbose=true, window=7, shift=-0.)
    # errors78_4 = get_gps_error(p78*"velocities_2018_current.h5", all_stations, dset="velos_shifted/1", verbose=true, window=7, shift=-0.)

    for f in ["velocities_2018_stackavg_max800.h5", "velocities_2018_current.h5", "velocities_2018_stackavg_max800.h5", "velocities_2017_current.h5"]
        push!(errors, get_gps_error(p85*f, all_stations, dset="velos_shifted/1", verbose=true, window=7, shift=-0.))
        push!(errors, get_gps_error(p78*f, all_stations, dset="velos_shifted/1", verbose=true, window=7, shift=-0.))
    end
    df = DataFrame(all_stations, round.(errors, digits=2)...)
    show(stdout, MIME("text/latex"), df)
end
