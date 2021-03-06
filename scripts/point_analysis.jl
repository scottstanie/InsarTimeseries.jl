import PyPlot
plt = PyPlot
include("./testl1.jl")
include("./plotting.jl")

p2c = InsarTimeseries.PHASE_TO_CM
p2mm = InsarTimeseries.PHASE_TO_CM * 365 * 10

load_geolist_intlist = InsarTimeseries.load_geolist_intlist
prune_cor = InsarTimeseries.prune_cor
remove_outliers = InsarTimeseries.remove_outliers
shrink_baseline = InsarTimeseries.shrink_baseline

invert(B::AbstractArray{T,2}, u::AbstractArray{T,1}) where {T<:Number} =
    InsarTimeseries.invert_pixel_l1(u, B)
# invt(B, v) = [p2mm * (B \ v) ; p2mm * invert(B, v) ]
# gives the stacked version, sbas version, and L1 version
invt(B, v) = [p2mm * (sum(v) / sum(B));
              p2mm * (B \ v);
              p2mm * invert(B, v) ]

function prune_igrams(g, i, u, ns = 3, ct = nothing, cp = nothing, shrink = false)
    g, i, u = remove_outliers(g, i, u, mean_sigma_cutoff = ns)
    # i, u = prune_cor(i, u, cor_pixel=cp, cor_thresh=ct)
    i, u = shrink ? shrink_baseline(g, i, u) : (i, u)
    return g, i, u
end
function prunesolve(
    g,
    i,
    u,
    B,
    nsigma = 3;
    cor_thresh = nothing,
    cor_pixel = nothing,
    shrink = false,
)
    constant = true
    geo, ints, unws = prune_igrams(g, i, u, nsigma, cor_thresh, cor_pixel, shrink)
    println("Num geos removed: $(length(g) - length(geo))")
    B = InsarTimeseries.prepB(geo, ints, constant)
    invt(B, unws)
end


function demo_point(
    rowcol;
    sigma = 4,
    max_temp = 800,
    show = true,
    refpoint = nothing,
    max_date = nothing,
    unwfile = "unw_stack.h5",
    input_dset = "stack_flat_shifted",
    ignore_file = "geolist_ignore.txt",
)
    geolist, intlist, igram_idxs =
        load_geolist_intlist(unwfile, ignore_file, max_temp, max_date = max_date)
    B = InsarTimeseries.build_B_matrix(geolist, intlist)
    Blin = sum(B, dims = 2)
    geolist, intlistall, igram_idxs_all =
        load_geolist_intlist(unwfile, ignore_file, 2000, max_date = max_date)
    Ball = InsarTimeseries.build_B_matrix(geolist, intlistall)
    Blinall = sum(Ball, dims = 2)

    unw_vals = get_stack_vals(unwfile, rowcol..., 1, input_dset, igram_idxs)
    # cc_vals = get_stack_vals("cc_stack.h5", rowcol..., 1, "stack", igram_idxs);

    unw_valsall = get_stack_vals(unwfile, rowcol..., 1, input_dset, igram_idxs_all)
    # cc_valsall = get_stack_vals("cc_stack.h5", rowcol..., 1, "stack", igram_idxs_all);

    if !isnothing(refpoint)
        println("re shifting to use $refpoint as reference")
        refunw = get_stack_vals(unwfile, refpoint..., 1, input_dset, igram_idxs)
        refunwall = get_stack_vals(unwfile, refpoint..., 1, input_dset, igram_idxs_all)
        unw_vals .-= refunw
        unw_valsall .-= refunwall
    end

    # println("Solving: Shrink=true, max temp = $max_temp")
    # @show prunesolve(geolist, intlist, unw_vals, Blin, sigma, shrink=true)
    println("Solving: no outlier remove, max temp = $max_temp")
    @show prunesolve(geolist, intlist, unw_vals, Blin, 1000, shrink = false)

    println("Solving: Shrink=false, max temp = $max_temp")
    @show prunesolve(geolist, intlist, unw_vals, Blin, sigma, shrink = false)

    println("Solving: Shrink=false, all igrams")
    @show prunesolve(geolist, intlistall, unw_valsall, Blinall, sigma, shrink = false)
    !show && return

    # Scatter plot of unw vals vs baseline
    plt.figure()
    plt.scatter(Blinall, p2c * unw_valsall, label = "all")
    plt.scatter(Blin, p2c * unw_vals, label = "600 days or less")
    plt.xlabel("Baseline (days")
    plt.ylabel("CM")
    plt.title("unw vs baseline")

    # plt.figure()
    # plt.scatter(Blinall, cc_valsall, label="all")
    # plt.scatter(Blin, cc_vals, label="$max_temp days or less")
    # plt.xlabel("Baseline (days")
    # plt.ylabel("CM")
    # plt.title("correlation vs baseline")

    # Plot unregularized
    plt.figure()
    unreg =
        p2c .* InsarTimeseries.integrate_velocities(
            B \ unw_vals,
            InsarTimeseries.day_diffs(geolist),
        )
    plt.plot(geolist, unreg, label = "unreg", "bo-")
    plt.title("Unregularized solution")
    plt.ylabel("CM")

    plot_grouped_by_day(geolist, intlist, unw_vals, sigma)
    plot_big_days(geolist, intlist, unw_vals, nsigma = sigma)
    return unreg
end
