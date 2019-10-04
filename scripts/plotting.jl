using HDF5
import PyPlot
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

