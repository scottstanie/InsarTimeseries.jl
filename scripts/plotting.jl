using HDF5
import PyPlot
import Polynomials
plt = PyPlot

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


function plot_big_days(geolist, intlist, vals, B; to_cm=true, label=nothing, nsigma=1, color=nothing)

    # means = InsarTimeseries.mean_abs_val(geolist, intlist, v);
    # high_days = sort(collect(zip(means, geolist)), rev=true)[1:5]
    high_days = InsarTimeseries.nsigma_days(geolist, intlist, vals, nsigma)
    return plot_days(geolist, intlist, vals, B, high_days, to_cm=to_cm, label=label, color=color)
end

function plot_days(geolist, intlist, vals, B, date_arr; to_cm=true, label=nothing, color=nothing)
    scale = to_cm ? InsarTimeseries.PHASE_TO_CM : 1
    v = vals .* scale   

    plt.figure()
    plt.scatter(B, v, label=label)

    ym = maximum(abs.(v)) * 1.2
    ylims = to_cm ? (-ym, ym) : (0, ym)  # correlation: 0 as bottom

    # for ii in 1:3
    for (ii, dd) in enumerate(date_arr)
        # hh = highest_days[ii][2]
        plt.scatter(InsarTimeseries.Blins_by_date(dd, intlist, B),
                                         InsarTimeseries.vals_by_date(dd, intlist, v),
                                         color=color,
                                         label="date=$dd")
        plt.ylim(ylims)
    end
    plt.legend()
    plt.show(block=false)
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

    geolist, intlist500, valid_igram_indices500 = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 500)

    m, n = 3, 5
    fig, axes = plt.subplots(m, n)
    for ii in 1:m
        for jj in 1:n
            idx = (ii-1)*n + jj
            if idx > 14
                continue
            end
            name = station_name_list[idx]

        unw_vals = get_stack_vals("unw_stack.h5", name, 1, "stack_flat_shifted", valid_igram_indices500)

        g2, i2, u2 = InsarTimeseries.remove_outliers(geolist, intlist500, unw_vals)

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
