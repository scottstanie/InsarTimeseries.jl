import PyPlot
plt = PyPlot
include("./testl1.jl")

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
    geolist, intlist500, valid_igram_indices500 = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 500)
    B500 = InsarTimeseries.build_B_matrix(geolist, intlist500);
    Blin500 = sum(B500, dims=2);
    unw_vals500_subs = get_stack_vals("unw_stack.h5", rowcol..., 1, "stack_flat_shifted", valid_igram_indices500)
    cc500_subs = get_stack_vals("cc_stack.h5", rowcol..., 1, "stack", valid_igram_indices500);
    prunesolve(geolist, intlistall, unw_valsall_subs, Blinall, 3, shrink=true)
    prunesolve(geolist, intlistall, unw_valsall_subs, Blinall, 3, shrink=false)
end
