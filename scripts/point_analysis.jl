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

invert(B::AbstractArray{T, 2}, u::AbstractArray{T, 1}) where {T<: Number} = InsarTimeseries.invert_pixel_l1(u, B)
invt(B, v) = [p2mm * (B \ v) ; p2mm * invert(B, v) ]

function prune_igrams(g, i, u, ns=3, ct=nothing, cp=nothing, shrink=true)
    g, i, u = remove_outliers(g, i, u, mean_sigma_cutoff=ns)
    # i, u = prune_cor(i, u, cor_pixel=cp, cor_thresh=ct)
    i, u = shrink ? shrink_baseline(g, i, u) : (i, u)
    return g, i, u
end
function prunesolve(g, i, u, B, nsigma=3; cor_thresh=nothing, cor_pixel=nothing, shrink=true)
    constant = true
    geo, ints, unws = prune_igrams(g, i, u, nsigma, cor_thresh, cor_pixel, shrink)
    B = InsarTimeseries.prepB(geo, ints, constant)
    invt(B, unws)
end


function demo_point(rowcol)
    geolist, intlist600, valid_igram_indices600 = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 600)
    B600 = InsarTimeseries.build_B_matrix(geolist, intlist600);
    Blin600 = sum(B600, dims=2);
    geolist, intlistall, valid_igram_indicesall = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 2000)
    Ball = InsarTimeseries.build_B_matrix(geolist, intlistall);
    Blinall = sum(Ball, dims=2);

    unw_vals600 = get_stack_vals("unw_stack.h5", rowcol..., 1, "stack_flat_shifted", valid_igram_indices600)
    cc_vals600 = get_stack_vals("cc_stack.h5", rowcol..., 1, "stack", valid_igram_indices600);

    unw_valsall = get_stack_vals("unw_stack.h5", rowcol..., 1, "stack_flat_shifted", valid_igram_indicesall)
    cc_valsall = get_stack_vals("cc_stack.h5", rowcol..., 1, "stack", valid_igram_indicesall);

    prunesolve(geolist, intlist600, unw_vals600, Blin600, 3, shrink=true)
    prunesolve(geolist, intlistall, unw_valsall, Blinall, 3, shrink=false)

    # Scatter plot of unw vals vs baseline
     plt.figure()
     plt.scatter(Blinall, p2c * unw_valsall, label="all")
     plt.scatter(Blin600, p2c * unw_vals600, label="600 days or less")
     plt.xlabel("Baseline (days")
     plt.ylabel("CM")
     plt.title("unw vs baseline")

     plt.figure()
     plt.scatter(Blinall, cc_valsall, label="all")
     plt.scatter(Blin600, cc_vals600, label="600 days or less")
     plt.xlabel("Baseline (days")
     plt.ylabel("CM")
     plt.title("correlation vs baseline")

     # Plot unregularized
     plt.figure()
     plt.plot(geolist, p2c .* InsarTimeseries.integrate_velocities(B600 \ unw_vals600, InsarTimeseries.day_diffs(geolist)), label="unreg")
     plt.title("Unregularized solution")
     plt.ylabel("CM")

     plot_grouped_by_day(geolist, intlist600, unw_vals600)
     plot_big_days(geolist, intlist600, unw_vals600, Blin600, nsigma=3)
end
