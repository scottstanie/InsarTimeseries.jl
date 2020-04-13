import PyPlot
plt = PyPlot
import InsarTimeseries
include("./testl1.jl")
include("./point_analysis.jl")


rowcol2 = [245, 189]  # small/looked uplift
# rowcol2 = [1224, 944]  # full uplift
# rowcol3 = [368, 131]  # Small/looked subs
# rowcol3 = [2229, 2236]  # TESTING for diffs
rowcol3 = [312, 135]  # Subs line that disappears with long baseline
# rowcol3 = [920, 329]  # Subs line that disappears with long baseline (big map)
rowcol4 = [32, 352]  # Subs around the decorrelated top

# rowcol = [1632, 687]  # Test near pecos, disappears with too much filtering
# rowcol2 = [1478, 743]  # bowl of strong subsidence that appears in 400-day soln, not in 600+
# rowcol3 = [1649, 699]  # patch of subsidence that appears in 400-day soln, not in 600+

# geolistig, intlist500ig, valid_igram_indices500ig = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 500)
geolist, intlist400, valid_igram_indices400 =
    load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 400)
_, intlist500, valid_igram_indices500 =
    load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 500)
_, intlist600, valid_igram_indices600 =
    load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 600)
_, intlist700, valid_igram_indices700 =
    load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 700)
_, intlist800, valid_igram_indices800 =
    load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 800)
_, intlistall, valid_igram_indicesall =
    load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 1500)
# geolist, intlist400, valid_igram_indices400 = load_geolist_intlist("unw_stack.h5", nothing, 400)


B400 = InsarTimeseries.build_B_matrix(geolist, intlist400);
Blin400 = sum(B400, dims = 2);
B500 = InsarTimeseries.build_B_matrix(geolist, intlist500);
Blin500 = sum(B500, dims = 2);
Blin600 = sum(InsarTimeseries.build_B_matrix(geolist, intlist600), dims = 2);
Blin700 = sum(InsarTimeseries.build_B_matrix(geolist, intlist700), dims = 2);
Blin800 = sum(InsarTimeseries.build_B_matrix(geolist, intlist800), dims = 2);
Blinall = sum(InsarTimeseries.build_B_matrix(geolist, intlistall), dims = 2);

println("point 0")

unw_vals500_uplift = get_stack_vals(
    "unw_stack.h5",
    rowcol2...,
    1,
    "stack_flat_shifted",
    valid_igram_indices500,
    reference_station = nothing,
);
unw_valsall_uplift = get_stack_vals(
    "unw_stack.h5",
    rowcol2...,
    1,
    "stack_flat_shifted",
    valid_igram_indicesall,
);
cc500_uplift = get_stack_vals("cc_stack.h5", rowcol2..., 1, "stack", valid_igram_indices500);

println("point 1")

unw_vals500_subs = get_stack_vals(
    "unw_stack.h5",
    rowcol3...,
    1,
    "stack_flat_shifted",
    valid_igram_indices500,
)
cc500_subs = get_stack_vals("cc_stack.h5", rowcol3..., 1, "stack", valid_igram_indices500);
unw_valsall_subs = get_stack_vals(
    "unw_stack.h5",
    rowcol3...,
    1,
    "stack_flat_shifted",
    valid_igram_indicesall,
);
ccall_subs = get_stack_vals("cc_stack.h5", rowcol3..., 1, "stack", valid_igram_indicesall);

println("point 2")
unw_valsall_txoz =
    get_stack_vals("unw_stack.h5", "TXOZ", 1, "stack_flat_shifted", valid_igram_indicesall);
unw_vals500_txoz = get_stack_vals(
    "unw_stack.h5",
    "TXOZ",
    1,
    "stack_flat_shifted",
    valid_igram_indices500,
    reference_station = nothing,
);
cc500_txoz = get_stack_vals("cc_stack.h5", "TXOZ", 1, "stack", valid_igram_indices500);
unw_vals500_nmhb = get_stack_vals(
    "unw_stack.h5",
    "NMHB",
    1,
    "stack_flat_shifted",
    valid_igram_indices500,
    reference_station = nothing,
);
cc500_nmhb = get_stack_vals("cc_stack.h5", "NMHB", 1, "stack", valid_igram_indices500);

println("point 3")
unw_vals500_txmh = get_stack_vals(
    "unw_stack.h5",
    "TXMH",
    1,
    "stack_flat_shifted",
    valid_igram_indices500,
    reference_station = nothing,
);
cc500_txmh = get_stack_vals("cc_stack.h5", "TXMH", 1, "stack", valid_igram_indices500);
cc500_txkm = get_stack_vals("cc_stack.h5", "TXKM", 1, "stack", valid_igram_indices500);
# 
println("done ")

_, intlist600, valid_igram_indices600 =
    load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 600)
B600 = InsarTimeseries.build_B_matrix(geolist, intlist600);
Blin600 = sum(B600, dims = 2);
unw_vals600_uplift = get_stack_vals(
    "unw_stack.h5",
    rowcol2...,
    1,
    "stack_flat_shifted",
    valid_igram_indices600,
    reference_station = nothing,
);
cc600_uplift = get_stack_vals(
    "cc_stack.h5",
    rowcol2...,
    1,
    "stack",
    valid_igram_indices600,
    reference_station = nothing,
);
unw_vals600_subs = get_stack_vals(
    "unw_stack.h5",
    rowcol3...,
    1,
    "stack_flat_shifted",
    valid_igram_indices600,
    reference_station = nothing,
);
cc600_subs = get_stack_vals(
    "cc_stack.h5",
    rowcol3...,
    1,
    "stack",
    valid_igram_indices600,
    reference_station = nothing,
);
unw_vals600_txoz = get_stack_vals(
    "unw_stack.h5",
    "TXOZ",
    1,
    "stack_flat_shifted",
    valid_igram_indices600,
    reference_station = nothing,
);
cc600_txoz = get_stack_vals(
    "cc_stack.h5",
    "TXOZ",
    1,
    "stack",
    valid_igram_indices600,
    reference_station = nothing,
);
unw_vals600_nmhb = get_stack_vals(
    "unw_stack.h5",
    "NMHB",
    1,
    "stack_flat_shifted",
    valid_igram_indices600,
    reference_station = nothing,
);
cc600_nmhb = get_stack_vals(
    "cc_stack.h5",
    "NMHB",
    1,
    "stack",
    valid_igram_indices600,
    reference_station = nothing,
);

function get_vals(days; station = nothing, rowcol = nothing, ignore = true)
    gi = ignore ? "geolist_ignore.txt" : nothing

    gl, il, vv = load_geolist_intlist("unw_stack.h5", gi, days)
    Bl = sum(InsarTimeseries.build_B_matrix(gl, il), dims = 2)
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
# rowcol_top = [25, 47]  # Uplift bowl in top left
rowcol_top = [14, 66]  # Uplift Edge of top left (random point)

unw_vals500_asc = get_stack_vals(
    "unw_stack.h5",
    rowcol_asc...,
    1,
    "stack_flat_shifted",
    valid_igram_indices500,
    reference_station = nothing,
);
cc500_asc = get_stack_vals(
    "cc_stack.h5",
    rowcol_asc...,
    1,
    "stack",
    valid_igram_indices500,
    reference_station = nothing,
);
unw_vals600_asc = get_stack_vals(
    "unw_stack.h5",
    rowcol_asc...,
    1,
    "stack_flat_shifted",
    valid_igram_indices600,
    reference_station = nothing,
);
# unw_vals500_desc = get_stack_vals("unw_stack.h5", rowcol_desc..., 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);

unw_vals500_top = get_stack_vals(
    "unw_stack.h5",
    rowcol_top...,
    1,
    "stack_flat_shifted",
    valid_igram_indices500,
);
cc500_top = get_stack_vals("cc_stack.h5", rowcol_top..., 1, "stack", valid_igram_indices500);

# g35, i35, v35 = load_geolist_intlist("unw_stack.h5", "geolist_ignore_35pct.txt", 500)
# Blin35 = sum(InsarTimeseries.build_B_matrix(g35, i35), dims=2)
# g50, i50, v50 = load_geolist_intlist("unw_stack.h5", "geolist_ignore_50pct.txt", 500)
# Blin50 = sum(InsarTimeseries.build_B_matrix(g50, i50), dims=2)
# unw_vals_up1 = get_stack_vals("unw_stack.h5", rowcol..., 1, "stack_flat_shifted", valid_igram_indices500);
# unw_vals_up35 = get_stack_vals("unw_stack.h5", rowcol..., 1, "stack_flat_shifted", v35);
# unw_vals_up50 = get_stack_vals("unw_stack.h5", rowcol..., 1, "stack_flat_shifted", v50);



#### 
# # Splitting by 2017 stuff
# geolist1, intlist1, valid_igram_indices1 = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 500, max_date=Date(2017,1,1))
# geolist2, intlist2, valid_igram_indices2 = load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 500, min_date=Date(2017,1,1))
# 
# B1 = InsarTimeseries.build_B_matrix(geolist1, intlist1)
# Blin1 = sum(B1, dims=2);
# B2 = InsarTimeseries.build_B_matrix(geolist2, intlist2)
# Blin2 = sum(B2, dims=2);
# 
# unw_vals1_txoz = get_stack_vals("unw_stack.h5", "TXOZ", 1, "stack_flat_shifted", valid_igram_indices1, reference_station=nothing);
# cc1_txoz = get_stack_vals("cc_stack.h5", "TXOZ", 1, "stack", valid_igram_indices1, reference_station=nothing);
# 
# unw_vals2_txoz = get_stack_vals("unw_stack.h5", "TXOZ", 1, "stack_flat_shifted", valid_igram_indices2, reference_station=nothing);
# cc2_txoz = get_stack_vals("cc_stack.h5", "TXOZ", 1, "stack", valid_igram_indices2, reference_station=nothing);
# 
# # Good area: txmh
# unw_vals1_txmh = get_stack_vals("unw_stack.h5", "TXMH", 1, "stack_flat_shifted", valid_igram_indices1, reference_station=nothing);
# cc1_txmh = get_stack_vals("cc_stack.h5", "TXMH", 1, "stack", valid_igram_indices1, reference_station=nothing);
# 
# unw_vals2_txmh = get_stack_vals("unw_stack.h5", "TXMH", 1, "stack_flat_shifted", valid_igram_indices2, reference_station=nothing);
# cc2_txmh = get_stack_vals("cc_stack.h5", "TXMH", 1, "stack", valid_igram_indices2, reference_station=nothing);
# 
# # TXS3: worst for 2017 only
# unw_vals1_txs3 = get_stack_vals("unw_stack.h5", "TXS3", 1, "stack_flat_shifted", valid_igram_indices1, reference_station=nothing);
# cc1_txs3 = get_stack_vals("cc_stack.h5", "TXS3", 1, "stack", valid_igram_indices1, reference_station=nothing);
# 
# unw_vals2_txs3 = get_stack_vals("unw_stack.h5", "TXS3", 1, "stack_flat_shifted", valid_igram_indices2, reference_station=nothing);
# cc2_txs3 = get_stack_vals("cc_stack.h5", "TXS3", 1, "stack", valid_igram_indices2, reference_station=nothing);

function plot_unregs(; sigma = 3, max_temp = 700, shrink = false, max_date = nothing)
    constant = false
    g, i, igram_idxs = load_geolist_intlist(
        "unw_stack.h5",
        "geolist_ignore.txt",
        max_temp,
        max_date = max_date,
    )
    rowcols = [
        [(414, 678), "noise? -8/yr for sigma=4+, 0 on sigma 3"],
        [(410, 667), "No signal, lots of stuff removed"],
        [(363, 594), "Bot right (real) bowl"],
        [(301, 301), "middle of map, vague red"],
        [(257, 81), "Bowl with eastward near pecos"],
        [(365, 128), "Deepest part of pecos wing"],
    ]
    plt.figure()
    for (rowcol, label) in rowcols
        row, col = rowcol
        unw_vals =
            get_stack_vals("unw_stack.h5", row, col, 1, "stack_flat_shifted", igram_idxs)

        geo, ints, unws = prune_igrams(g, i, unw_vals, sigma, nothing, nothing, false)
        ints, unws = shrink ? shrink_baseline(geo, ints, unws) : (ints, unws)
        B = InsarTimeseries.prepB(geo, ints, constant)
        soln =
            p2c .*
            InsarTimeseries.integrate_velocities(B \ unws, InsarTimeseries.day_diffs(geo))

        Blin = sum(B, dims = 2)
        soln_lin = (Blin\unws)[1] * p2mm
        println("$label: $soln_lin cm/yr")

        plt.plot(
            geo,
            soln,
            ".-",
            label = "($row, $col): $label, Linear: $(round(soln_lin)) mm/yr",
        )
    end
    plt.legend()
    plt.title("Sigma = $sigma, max baseline = $max_temp, shrink = $shrink")
end
