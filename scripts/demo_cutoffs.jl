import PyPlot
plt = PyPlot
import InsarTimeseries
include("./testl1.jl")

invert(B::AbstractArray{T, 2}, u::AbstractArray{T, 1}) where {T<: Number} = InsarTimeseries.invert_pixel(u, B, rho=1.0, alpha=1.5)
invt(B, v) = [p2mm * (B \ v) ; p2mm * invert(B, v) ]

# rowcol2 = [245, 189]  # small/looked uplift
rowcol2 = [1224, 944]  # full uplift
# rowcol3 = [368, 131]  # Small/looked subs
rowcol3 = [2229, 2236]  # TESTING for diffs

rowcol = [1632, 687]  # Test near pecos, disappears with too much filtering
# rowcol2 = [1478, 743]  # bowl of strong subsidence that appears in 400-day soln, not in 600+
# rowcol3 = [1649, 699]  # patch of subsidence that appears in 400-day soln, not in 600+

# geolistig, intlist400ig, valid_igram_indices400ig = InsarTimeseries.load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 400)
geolist, intlist400, valid_igram_indices400 = InsarTimeseries.load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 400)
_, intlist600, valid_igram_indices600 = InsarTimeseries.load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 600)
# geolist, intlist400, valid_igram_indices400 = InsarTimeseries.load_geolist_intlist("unw_stack.h5", nothing, 400)
B400 = InsarTimeseries.build_B_matrix(geolist, intlist400);
Blin400 = sum(B400, dims=2);
Blin600 = sum(InsarTimeseries.build_B_matrix(geolist, intlist600), dims=2);

unw_vals400_uplift = get_stack_vals("unw_stack.h5", rowcol2..., 1, "stack_flat_shifted", valid_igram_indices400, reference_station=nothing);
cc400_uplift = get_stack_vals("cc_stack.h5", rowcol2..., 1, "stack", valid_igram_indices400, reference_station=nothing);
unw_vals400_subs = get_stack_vals("unw_stack.h5", rowcol..., 1, "stack_flat_shifted", valid_igram_indices400, reference_station=nothing);
cc400_subs = get_stack_vals("cc_stack.h5", rowcol..., 1, "stack", valid_igram_indices400, reference_station=nothing);
unw_vals400_txoz = get_stack_vals("unw_stack.h5", "TXOZ", 1, "stack_flat_shifted", valid_igram_indices400, reference_station=nothing);
cc400_txoz = get_stack_vals("cc_stack.h5", "TXOZ", 1, "stack", valid_igram_indices400, reference_station=nothing);
unw_vals400_nmhb = get_stack_vals("unw_stack.h5", "NMHB", 1, "stack_flat_shifted", valid_igram_indices400, reference_station=nothing);
cc400_nmhb = get_stack_vals("cc_stack.h5", "NMHB", 1, "stack", valid_igram_indices400, reference_station=nothing);

# 
#
# TODO: CHANGE BACK 700 to 600
_, intlist600, valid_igram_indices600 = InsarTimeseries.load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 600)
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

    gl, il, vv = InsarTimeseries.load_geolist_intlist("unw_stack.h5", gi, days)
    Bl = sum(InsarTimeseries.build_B_matrix(gl, il), dims=2);
    p = isnothing(station) ? rowcol : [station]
    println(p)

    u = get_stack_vals("unw_stack.h5", p..., 1, "stack_flat_shifted", vv)
    c = get_stack_vals("cc_stack.h5", p..., 1, "stack", vv)
    return gl, il, vv, Bl, u, c
end

#### Points for looked at the diffs:
# Looked, path78 asc
rowcol_asc = [226, 40]
# # Looked, path85 desc
# rowcol_desc = [106, 400]

unw_vals400_asc = get_stack_vals("unw_stack.h5", rowcol_asc..., 1, "stack_flat_shifted", valid_igram_indices400, reference_station=nothing);
cc400_asc = get_stack_vals("cc_stack.h5", rowcol_asc..., 1, "stack", valid_igram_indices400, reference_station=nothing);
unw_vals600_asc = get_stack_vals("unw_stack.h5", rowcol_asc..., 1, "stack_flat_shifted", valid_igram_indices600, reference_station=nothing);
# unw_vals400_desc = get_stack_vals("unw_stack.h5", rowcol_desc..., 1, "stack_flat_shifted", valid_igram_indices400, reference_station=nothing);

