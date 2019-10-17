import PyPlot
plt = PyPlot
import InsarTimeseries
include("./testl1.jl")

load_geolist_intlist = InsarTimeseries.load_geolist_intlist
prune_cor = InsarTimeseries.prune_cor
remove_outliers = InsarTimeseries.remove_outliers
shrink_baseline = InsarTimeseries.shrink_baseline

invert(B::AbstractArray{T, 2}, u::AbstractArray{T, 1}) where {T<: Number} = InsarTimeseries.invert_pixel_l1(u, B)
invt(B, v) = [p2mm * (B \ v) ; p2mm * invert(B, v) ]

function prune_igrams(g, i, u, ns=3, ct=nothing, cp=nothing)
    g, i, u = remove_outliers(g, i, u, mean_sigma_cutoff=ns)
    # i, u = prune_cor(i, u, cor_pixel=cp, cor_thresh=ct)
    i, u = shrink_baseline(g, i, u)
    return g, i, u
end

function prunesolve(g, i, u, B, nsigma=3, cor_thresh=nothing, cor_pixel=nothing)
    constant = true
    geo, ints, unws = prune_igrams(g, i, u, nsigma, cor_thresh, cor_pixel)
    B = InsarTimeseries.prepB(geo, ints, constant)
    invt(B, unws)
end


rowcol2 = [245, 189]  # small/looked uplift
# rowcol2 = [1224, 944]  # full uplift
rowcol3 = [368, 131]  # Small/looked subs
# rowcol3 = [2229, 2236]  # TESTING for diffs

# rowcol = [1632, 687]  # Test near pecos, disappears with too much filtering
# rowcol2 = [1478, 743]  # bowl of strong subsidence that appears in 400-day soln, not in 600+
# rowcol3 = [1649, 699]  # patch of subsidence that appears in 400-day soln, not in 600+

# geolistig, intlist500ig, valid_igram_indices500ig = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 500)
geolist, intlist400, valid_igram_indices400 = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 400)
_, intlist500, valid_igram_indices500 = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 500)
_, intlist600, valid_igram_indices600 = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 600)
# geolist, intlist400, valid_igram_indices400 = load_geolist_intlist("unw_stack.h5", nothing, 400)
B400 = InsarTimeseries.build_B_matrix(geolist, intlist400);
Blin400 = sum(B400, dims=2);
B500 = InsarTimeseries.build_B_matrix(geolist, intlist500);
Blin500 = sum(B500, dims=2);
Blin600 = sum(InsarTimeseries.build_B_matrix(geolist, intlist600), dims=2);

unw_vals500_uplift = get_stack_vals("unw_stack.h5", rowcol2..., 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);
cc500_uplift = get_stack_vals("cc_stack.h5", rowcol2..., 1, "stack", valid_igram_indices500, reference_station=nothing);
unw_vals500_subs = get_stack_vals("unw_stack.h5", rowcol3..., 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);
cc500_subs = get_stack_vals("cc_stack.h5", rowcol3..., 1, "stack", valid_igram_indices500, reference_station=nothing);
unw_vals500_txoz = get_stack_vals("unw_stack.h5", "TXOZ", 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);
cc500_txoz = get_stack_vals("cc_stack.h5", "TXOZ", 1, "stack", valid_igram_indices500, reference_station=nothing);
unw_vals500_nmhb = get_stack_vals("unw_stack.h5", "NMHB", 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);
cc500_nmhb = get_stack_vals("cc_stack.h5", "NMHB", 1, "stack", valid_igram_indices500, reference_station=nothing);

unw_vals500_txmh = get_stack_vals("unw_stack.h5", "TXMH", 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);
cc500_txmh = get_stack_vals("cc_stack.h5", "TXMH", 1, "stack", valid_igram_indices500, reference_station=nothing);
# 
#
_, intlist600, valid_igram_indices600 = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 600)
B600 = InsarTimeseries.build_B_matrix(geolist, intlist600);
Blin600 = sum(B600, dims=2);
unw_vals600_uplift = get_stack_vals("unw_stack.h5", rowcol2..., 1, "stack_flat_shifted", valid_igram_indices600, reference_station=nothing);
cc600_uplift = get_stack_vals("cc_stack.h5", rowcol2..., 1, "stack", valid_igram_indices600, reference_station=nothing);
unw_vals600_subs = get_stack_vals("unw_stack.h5", rowcol3..., 1, "stack_flat_shifted", valid_igram_indices600, reference_station=nothing);
cc600_subs = get_stack_vals("cc_stack.h5", rowcol3..., 1, "stack", valid_igram_indices600, reference_station=nothing);
unw_vals600_txoz = get_stack_vals("unw_stack.h5", "TXOZ", 1, "stack_flat_shifted", valid_igram_indices600, reference_station=nothing);
cc600_txoz = get_stack_vals("cc_stack.h5", "TXOZ", 1, "stack", valid_igram_indices600, reference_station=nothing);
unw_vals600_nmhb = get_stack_vals("unw_stack.h5", "NMHB", 1, "stack_flat_shifted", valid_igram_indices600, reference_station=nothing);
cc600_nmhb = get_stack_vals("cc_stack.h5", "NMHB", 1, "stack", valid_igram_indices600, reference_station=nothing);

function get_vals(days; station=nothing, rowcol=nothing, ignore=true)
    gi = ignore ? "geolist_ignore.txt" : nothing

    gl, il, vv = load_geolist_intlist("unw_stack.h5", gi, days)
    Bl = sum(InsarTimeseries.build_B_matrix(gl, il), dims=2);
    p = isnothing(station) ? rowcol : [station]
    println(p)

    u = get_stack_vals("unw_stack.h5", p..., 1, "stack_flat_shifted", vv)
    c = get_stack_vals("cc_stack.h5", p..., 1, "stack", vv)
    return gl, il, vv, Bl, u, c
end

#### Points for looked at the diffs:
# Looked, path78 asc
# rowcol_asc = [226, 40]
rowcol_asc = [259, 82]  # Uplift looks stronger in stack than not
# # Looked, path85 desc
# rowcol_desc = [106, 400]
#
#
# rowcol = [25, 47]  # Uplift bowl in top left
rowcol = [14, 66]  # Uplift Edge of top left (random point)

unw_vals500_asc = get_stack_vals("unw_stack.h5", rowcol_asc..., 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);
cc500_asc = get_stack_vals("cc_stack.h5", rowcol_asc..., 1, "stack", valid_igram_indices500, reference_station=nothing);
unw_vals600_asc = get_stack_vals("unw_stack.h5", rowcol_asc..., 1, "stack_flat_shifted", valid_igram_indices600, reference_station=nothing);
# unw_vals500_desc = get_stack_vals("unw_stack.h5", rowcol_desc..., 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);

g35, i35, v35 = load_geolist_intlist("unw_stack.h5", "geolist_ignore_35pct.txt", 500)
Blin35 = sum(InsarTimeseries.build_B_matrix(g35, i35), dims=2)
g50, i50, v50 = load_geolist_intlist("unw_stack.h5", "geolist_ignore_50pct.txt", 500)
Blin50 = sum(InsarTimeseries.build_B_matrix(g50, i50), dims=2)
unw_vals_up1 = get_stack_vals("unw_stack.h5", rowcol..., 1, "stack_flat_shifted", valid_igram_indices500);
unw_vals_up35 = get_stack_vals("unw_stack.h5", rowcol..., 1, "stack_flat_shifted", v35);
unw_vals_up50 = get_stack_vals("unw_stack.h5", rowcol..., 1, "stack_flat_shifted", v50);



#### 
# Splitting by 2017 stuff
geolist1, intlist1, valid_igram_indices1 = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 500, max_date=Date(2017,1,1))
geolist2, intlist2, valid_igram_indices2 = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 500, min_date=Date(2017,1,1))

B1 = InsarTimeseries.build_B_matrix(geolist1, intlist1)
Blin1 = sum(B1, dims=2);
B2 = InsarTimeseries.build_B_matrix(geolist2, intlist2)
Blin2 = sum(B2, dims=2);

unw_vals1_txoz = get_stack_vals("unw_stack.h5", "TXOZ", 1, "stack_flat_shifted", valid_igram_indices1, reference_station=nothing);
cc1_txoz = get_stack_vals("cc_stack.h5", "TXOZ", 1, "stack", valid_igram_indices1, reference_station=nothing);

unw_vals2_txoz = get_stack_vals("unw_stack.h5", "TXOZ", 1, "stack_flat_shifted", valid_igram_indices2, reference_station=nothing);
cc2_txoz = get_stack_vals("cc_stack.h5", "TXOZ", 1, "stack", valid_igram_indices2, reference_station=nothing);

# Good area: txmh
unw_vals1_txmh = get_stack_vals("unw_stack.h5", "TXMH", 1, "stack_flat_shifted", valid_igram_indices1, reference_station=nothing);
cc1_txmh = get_stack_vals("cc_stack.h5", "TXMH", 1, "stack", valid_igram_indices1, reference_station=nothing);

unw_vals2_txmh = get_stack_vals("unw_stack.h5", "TXMH", 1, "stack_flat_shifted", valid_igram_indices2, reference_station=nothing);
cc2_txmh = get_stack_vals("cc_stack.h5", "TXMH", 1, "stack", valid_igram_indices2, reference_station=nothing);

# TXS3: worst for 2017 only
unw_vals1_txs3 = get_stack_vals("unw_stack.h5", "TXS3", 1, "stack_flat_shifted", valid_igram_indices1, reference_station=nothing);
cc1_txs3 = get_stack_vals("cc_stack.h5", "TXS3", 1, "stack", valid_igram_indices1, reference_station=nothing);

unw_vals2_txs3 = get_stack_vals("unw_stack.h5", "TXS3", 1, "stack_flat_shifted", valid_igram_indices2, reference_station=nothing);
cc2_txs3 = get_stack_vals("cc_stack.h5", "TXS3", 1, "stack", valid_igram_indices2, reference_station=nothing);
