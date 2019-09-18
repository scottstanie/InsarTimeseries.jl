import PyPlot
plt = PyPlot
import InsarTimeseries
include("./testl1.jl")

invert(B::AbstractArray{T, 2}, u::AbstractArray{T, 1}) where {T<: Number} = InsarTimeseries.invert_pixel(u, B, rho=1.0, alpha=1.5)
invt(B, v) = [p2mm * (B \ v) ; p2mm * invert(B, v) ]

rowcol2 = [245, 189]  # small/looked uplift
# rowcol2 = [1224, 944]  # full uplift
rowcol3 = [368, 131]  # Small/looked subs
# rowcol3 = [2229, 2236]  # TESTING for diffs

# geolistig, intlist500ig, valid_igram_indices500ig = InsarTimeseries.load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 500)
geolist, intlist500, valid_igram_indices500 = InsarTimeseries.load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 500)
_, intlist300, valid_igram_indices300 = InsarTimeseries.load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 300)
# geolist, intlist500, valid_igram_indices500 = InsarTimeseries.load_geolist_intlist("unw_stack.h5", nothing, 500)
B500 = InsarTimeseries.build_B_matrix(geolist, intlist500);
Blin500 = sum(B500, dims=2);
Blin300 = sum(InsarTimeseries.build_B_matrix(geolist, intlist300), dims=2);

unw_vals500_uplift = get_stack_vals("unw_stack.h5", rowcol2..., 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);
cc500_uplift = get_stack_vals("cc_stack.h5", rowcol2..., 1, "stack", valid_igram_indices500, reference_station=nothing);
unw_vals500_subs = get_stack_vals("unw_stack.h5", rowcol3..., 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);
cc500_subs = get_stack_vals("cc_stack.h5", rowcol3..., 1, "stack", valid_igram_indices500, reference_station=nothing);
unw_vals500_txoz = get_stack_vals("unw_stack.h5", "TXOZ", 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);
cc500_txoz = get_stack_vals("cc_stack.h5", "TXOZ", 1, "stack", valid_igram_indices500, reference_station=nothing);
unw_vals500_nmhb = get_stack_vals("unw_stack.h5", "NMHB", 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);
cc500_nmhb = get_stack_vals("cc_stack.h5", "NMHB", 1, "stack", valid_igram_indices500, reference_station=nothing);

# 
#
# TODO: CHANGE BACK 700 to 300
_, intlist300, valid_igram_indices300 = InsarTimeseries.load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 700)
B300 = InsarTimeseries.build_B_matrix(geolist, intlist300);
Blin300 = sum(B300, dims=2);
unw_vals300_uplift = get_stack_vals("unw_stack.h5", rowcol2..., 1, "stack_flat_shifted", valid_igram_indices300, reference_station=nothing);
cc300_uplift = get_stack_vals("cc_stack.h5", rowcol2..., 1, "stack", valid_igram_indices300, reference_station=nothing);
unw_vals300_subs = get_stack_vals("unw_stack.h5", rowcol3..., 1, "stack_flat_shifted", valid_igram_indices300, reference_station=nothing);
cc300_subs = get_stack_vals("cc_stack.h5", rowcol3..., 1, "stack", valid_igram_indices300, reference_station=nothing);
unw_vals300_txoz = get_stack_vals("unw_stack.h5", "TXOZ", 1, "stack_flat_shifted", valid_igram_indices300, reference_station=nothing);
cc300_txoz = get_stack_vals("cc_stack.h5", "TXOZ", 1, "stack", valid_igram_indices300, reference_station=nothing);
unw_vals300_nmhb = get_stack_vals("unw_stack.h5", "NMHB", 1, "stack_flat_shifted", valid_igram_indices300, reference_station=nothing);
cc300_nmhb = get_stack_vals("cc_stack.h5", "NMHB", 1, "stack", valid_igram_indices300, reference_station=nothing);

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

unw_vals500_asc = get_stack_vals("unw_stack.h5", rowcol_asc..., 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);
cc500_asc = get_stack_vals("cc_stack.h5", rowcol_asc..., 1, "stack", valid_igram_indices500, reference_station=nothing);
unw_vals300_asc = get_stack_vals("unw_stack.h5", rowcol_asc..., 1, "stack_flat_shifted", valid_igram_indices300, reference_station=nothing);
# unw_vals500_desc = get_stack_vals("unw_stack.h5", rowcol_desc..., 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);

