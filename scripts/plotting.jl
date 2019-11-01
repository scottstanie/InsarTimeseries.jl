using HDF5
import PyPlot
import Polynomials
import Dierckx: Spline2D
plt = PyPlot
include("./colors.jl")
 
load_geolist_intlist = InsarTimeseries.load_geolist_intlist
baseline = InsarTimeseries.temporal_baseline

function save_paper_figure(fig, fname, axis_off=false)
    fig.tight_layout()
    axis_off && plt.axis("off")
    println("Saving $fname")
    fig.savefig(fname, bbox_inches="tight", transparent=true, dpi=300)
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

    fig, axes = plt.subplots(1, n, squeeze=false)
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
    end
    fig.legend()
end

function plot_big_days(geolist, intlist, vals, B; to_cm=true, label=nothing, nsigma=3, min_spread=0, color=nothing)

    # means = InsarTimeseries.mean_abs_val(geolist, intlist, v);
    # high_days = sort(collect(zip(means, geolist)), rev=true)[1:5]
    high_days = InsarTimeseries.nsigma_days(geolist, intlist, vals, nsigma, min_spread)
    return plot_days(geolist, intlist, vals, B, high_days, to_cm=to_cm, label=label, color=color)
end

function plot_days(geolist, intlist, vals, B, date_arr; to_cm=true, label=nothing, color=nothing)
    scale = to_cm ? InsarTimeseries.PHASE_TO_CM : 1
    v = vals .* scale   

    fig, ax =plt.subplots()
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
    fig.legend()
    plt.show(block=false)
    return fig, ax
end


# Plot gps east up data
function plot_eu(station_name, insar_slopes=nothing, marker=".")
    dts, east, north, up = InsarTimeseries.get_gps_enu(station_name)

    insar_east, insar_up = !isnothing(insar_slopes) ? insar_slopes : (nothing, nothing)

    fig, axes = plt.subplots(1, 2)
    axes[1].plot(dts, east, marker, label="gps")
    axes[1].set_title("east")
    axes[1].set_ylabel("cm")
    axes[1].set_ylim((-2.5, 2.5))
    if !isnothing(insar_east)
        plot_insar(axes[1], dts, insar_east)
    end

    axes[2].plot(dts, up, marker, label="gps")
    axes[2].set_title("up")
    axes[2].set_ylim((-2.5, 2.5))
    if !isnothing(insar_up)
        plot_insar(axes[2], dts, insar_up)
    end

    fig.suptitle(station_name)
    return fig, axes
end

function plot_insar(ax, dts, insar_slope; label="insar")
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
                                   vm=15, bins=70)
    # geolist, excls = count_dates(fname, dset)
    # @time excl_map = InsarTimeseries.decode_excluded(excls, geolist);
    p2c = InsarTimeseries.PHASE_TO_CM

    geolist, intlist600, valid_igram_indices600 = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 600)

    m, n = 3, 5
    fig, axes = plt.subplots(m, n)
    for ii in 1:m
        for jj in 1:n
            idx = (ii-1)*n + jj
            if idx > 14
                continue
            end
            name = station_name_list[idx]

        unw_vals = get_stack_vals("unw_stack.h5", name, 1, "stack_flat_shifted", valid_igram_indices600)

        g2, i2, u2 = InsarTimeseries.remove_outliers(geolist, intlist600, unw_vals)

        ax = axes[ii, jj]

        n1, _, _ = ax.hist(unw_vals .* p2c, range=(-vm, vm), bins=bins, label="Original spread")
        n1, _, _ = ax.hist(u2 .* p2c, range=(-vm, vm), bins=bins, label="Outliers removed")
        ax.set_title("$name")
        ax.grid("on")
    end
    end
    plt.legend()
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
function plot_image_line(img::MapImages.MapImage, rowcol1, rowcol2; plotkwargs...)
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

function plot_gps_station(name, insar_mm; ref="TXKM", ylim=(-3, 3), title="")
    fig, ax = plt.subplots()
    dts, los = get_gps_los(name, reference_station=ref)
    day_nums = _get_day_nums(dts)

    ax.plot(dts, los, "b.", markersize=3, label="gps")

    insar_cm_day = insar_mm / 365 / 10

    full_defo = insar_cm_day * (dts[end] - dts[1]).value
    ax.plot(dts, -full_defo/2 .+ day_nums .* insar_cm_day, "r", lw=3)

    # Or if you want different settings for the grids:
    ax.grid(which="major", alpha=0.5)
    ax.set_xticks(dts[1]:Dates.Day(365):dts[end])
    ax.set_yticks([-2, 0, 2])

    ax.set_ylim(ylim)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 20
    rcParams["font.weight"] = "bold"
#     for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#              ax.get_xticklabels() + ax.get_yticklabels()):
#     item.set_fontsize(20)
    return fig
end


