import PyPlot
plt = PyPlot
include("./testl1.jl")

invt(B, v) = [p2mm * (B \ v) ; p2mm * invert(B, v) ]
invert(B::AbstractArray{T, 2}, u::AbstractArray{T, 1}) where {T<: Number} = InsarTimeseries.invert_pixel(u, B, rho=1.0, alpha=1.5)

rowcol2 = [245, 189]
rowcol3 = [368, 131]

geolist, intlist500, valid_igram_indices500 = InsarTimeseries.load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 500)
B500 = InsarTimeseries.build_B_matrix(geolist, intlist500);
Blin500 = sum(B500, dims=2);

unw_vals500_uplift = get_stack_vals("unw_stack.h5", rowcol2..., 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);
cc500_uplift = get_stack_vals("cc_stack.h5", rowcol2..., 1, "stack", valid_igram_indices500, reference_station=nothing);
unw_vals500_subs = get_stack_vals("unw_stack.h5", rowcol3..., 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);
cc500_subs = get_stack_vals("cc_stack.h5", rowcol3..., 1, "stack", valid_igram_indices500, reference_station=nothing);
unw_vals500_txoz = get_stack_vals("unw_stack.h5", "TXOZ", 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);
cc500_txoz = get_stack_vals("cc_stack.h5", "TXOZ", 1, "stack", valid_igram_indices500, reference_station=nothing);
unw_vals500_nmhb = get_stack_vals("unw_stack.h5", "NMHB", 1, "stack_flat_shifted", valid_igram_indices500, reference_station=nothing);
cc500_nmhb = get_stack_vals("cc_stack.h5", "NMHB", 1, "stack", valid_igram_indices500, reference_station=nothing);

_, intlist300, valid_igram_indices300 = InsarTimeseries.load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", 300)
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
