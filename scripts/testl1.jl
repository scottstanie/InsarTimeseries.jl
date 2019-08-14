using InsarTimeseries
using Convex
using Dates
using ECOS
using HDF5
using Plots
using PyCall
# using SCS

latlon = pyimport("apertools.latlon")
sario = pyimport("apertools.sario")
gps = pyimport("apertools.gps")

dem_rsc = sario.load("dem.rsc")
station_name = "NMHB"
lon, lat = gps.station_lonlat(station_name)

row, col = latlon.nearest_pixel(dem_rsc, lon=lon, lat=lat)

row = convert(Int, row)
col = convert(Int, col)
println("Getting data at $station_name: row $row, col $col, lon $lon, lat $lat")

dts, los_gps_data = gps.load_gps_los_data("../", station_name, start_year=2015, end_year=2018)
dts = convert(Array{Dates.Date, 1}, dts);


defo_filename = "deformation_unreg_maxtemp400.h5"
defo_linear_filename = "deformation_linear_maxtemp400.h5"
println("Defo filename $defo_filename, defo linear: $defo_linear_filename")

unw_stack_file="unw_stack.h5"
ignore_geo_file = "geolist_ignore.txt"

max_temporal_baseline = 400
geolist, intlist, valid_igram_indices = InsarTimeseries.load_geolist_intlist(unw_stack_file, ignore_geo_file, max_temporal_baseline);

timediffs = InsarTimeseries.day_diffs(geolist)

B = InsarTimeseries.build_B_matrix(geolist, intlist)
Blin= sum(B, dims=2)

# The values already solved for:
ts_linear = h5read(defo_linear_filename, "stack", (col, row, :))[1, 1, :];
ts_unreg = h5read(defo_filename, "stack", (col, row, :))[1, 1, :];

dset = "stack_deramped_1_shifted"
println("Loading $row, $col from $dset")
@time unw_vals_all = h5read(unw_stack_file, dset, (col, row, :))[1, 1, :];
unw_vals = unw_vals_all[valid_igram_indices];

v_linear_lstsq = Blin \ unw_vals
ts_linear_lstsq = InsarTimeseries.PHASE_TO_CM .* InsarTimeseries.integrate_velocities(v_linear_lstsq, timediffs)
v_unreg_lstsq = B \ unw_vals
ts_unreg_lstsq = InsarTimeseries.PHASE_TO_CM .* InsarTimeseries.integrate_velocities(v_unreg_lstsq, timediffs);


# solver = SCSSolver()
solver = ECOSSolver(verbose=0)


v1linear = Variable(1)
prob_linear = minimize(norm(Blin*v1linear - unw_vals, 1))

v1unreg = Variable(size(B, 2))
prob_unreg = minimize(norm(B*v1unreg - unw_vals, 1))


solve!(prob_linear, solver)
solve!(prob_unreg, solver)

# To check that l2 norm min matches the backslash result
# v2linear_test = Variable(1)
# prob_linear_test2 = minimize(norm(Blin*v2linear_test - unw_vals, 2))
# solve!(prob_linear_test2, solver)

ts_linear_l1 = InsarTimeseries.PHASE_TO_CM .* InsarTimeseries.integrate_velocities(Float32.([v1linear.value]), timediffs);
ts_unreg_l1 = InsarTimeseries.PHASE_TO_CM .* InsarTimeseries.integrate_velocities(Float32.(reshape(v1unreg.value, :)), timediffs);


# p1 = plot(geolist, ts_linear)
# plot!(p1, geolist, ts_unreg)
p1 = plot(geolist, ts_linear_lstsq, title="L2 least squares solution")
plot!(p1, geolist, ts_unreg_lstsq)
plot!(p1, geolist, ts_unreg_lstsq)
plot!(p1, dts, los_gps_data, line=nothing, marker=:x, markersize=1)

p2 = plot(geolist, ts_linear_l1, title="L1 norm minimization")
plot!(p2, geolist, ts_unreg_l1)
plot!(p2, dts, los_gps_data, line=nothing, marker=:x, markersize=1)

plot(p1, p2)
