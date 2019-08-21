push!(LOAD_PATH,"/home/scott/repos/InsarTimeseries.jl/src/")
import InsarTimeseries
import InsarTimeseries: PHASE_TO_CM
using Convex
using Dates
using ECOS
using HDF5
using LinearAlgebra
using Statistics: mean
using Plots
using PyCall
import Polynomials
# using SCS

# TODO: USE REFERENCE STATION AND this
# #df = create_insar_gps_df(geo_path, defo_filename=defo_filename, reference_station=reference_station)

YLIMS = (-12, 12)
station_name_list = ["NMHB", "TXAD", "TXBG", "TXBL", "TXCE", "TXFS", "TXKM", 
                     "TXL2", "TXMC", "TXMH", "TXOE", "TXOZ", "TXS3", "TXSO"]

latlon = pyimport("apertools.latlon")
sario = pyimport("apertools.sario")
gps = pyimport("apertools.gps")

# const DEFO_FILENAME = "deformation_unreg_maxtemp400.h5"
# const DEFO_FILENAME_LINEAR = "deformation_linear_maxtemp400.h5"
# const DEFO_FILENAME_L1 = "deformation_linear_maxtemp400_l1.h5"
# println("Defo filename $DEFO_FILENAME, defo linear: $DEFO_FILENAME_LINEAR")

const UNW_STACK_FILE="unw_stack.h5"
const DSET = "stack_deramped_1_shifted"
const DSET_NOSHIFT = "stack_deramped_1"
const ignore_geo_file = "geolist_ignore.txt"
const max_temporal_baseline = 400

GEOLIST, INTLIST, VALID_IGRAM_INDICES = InsarTimeseries.load_geolist_intlist(UNW_STACK_FILE, ignore_geo_file, max_temporal_baseline);
timediffs = InsarTimeseries.day_diffs(GEOLIST)


# For checking / plotting specific interestion points:
# Sinkhole
row, col = 1838, 1367
# well uplift
row, col = 1223, 943


# Functions for error evaluation
rms(arr) = sqrt(mean(arr.^2))
total_abs_error(arr) = sum(abs.(arr))


function get_unw_vals(unw_stack_file, station_name)
    dem_rsc = sario.load("dem.rsc")
    lon, lat = gps.station_lonlat(station_name)
    row, col = map(x-> convert(Int, x), 
                   latlon.nearest_pixel(dem_rsc, lon=lon, lat=lat))

    println("Getting data at $station_name: row $row, col $col, lon $lon, lat $lat")
    @time unw_vals = _read_unw(unw_stack_file, row, col)
    return unw_vals
end

function get_unw_vals(unw_stack_file, row, col)
    @time unw_vals = _read_unw(unw_stack_file, row, col)
    return unw_vals
end

function _read_unw(unw_stack_file, row, col)
    println("Loading $row, $col from $DSET")
    # Note: swapping the col and row in h5read since julia is col-major
    unw_vals_all = h5read(unw_stack_file, DSET, (col, row, :))[1, 1, :]
    unw_vals = unw_vals_all[VALID_IGRAM_INDICES]
end

function get_gps_los(station_name, geo_path="../"; reference_station=nothing)
    dts, gps_los_data = gps.load_gps_los_data(geo_path, station_name, start_year=2015, end_year=2018)
    dts = convert(Array{Dates.Date, 1}, dts)
    return dts, gps_los_data
end

_get_day_nums(dts) = [( d - dts[1]).value for d in dts]

function fit_line(dts, data)
    day_nums = _get_day_nums(dts)
    p = Polynomials.polyfit(day_nums, data, 1)
    # p(day_nums[end])
    return p
end

function solve_insar_ts(station_name::String)
    unw_vals = get_unw_vals(UNW_STACK_FILE, station_name)
    return solve_insar_ts(unw_vals)
end

function solve_insar_ts(row::Int, col::Int)
    unw_vals = get_unw_vals(UNW_STACK_FILE, row, col)
    return solve_insar_ts(unw_vals)
end

function solve_insar_ts(unw_vals::Array{<:AbstractFloat, 1})

    B = InsarTimeseries.build_B_matrix(GEOLIST, INTLIST)
    Blin = sum(B, dims=2)

    v_linear_lstsq = Blin \ unw_vals
    v_unreg_lstsq = B \ unw_vals

    # # solver = SCSSolver()
    # solver = ECOSSolver(verbose=0)
    # prob_linear = minimize(norm(Blin*v_linear_l1 - unw_vals, 1))
    # prob_unreg = minimize(norm(B*v_unreg_l1 - unw_vals, 1))
    # solve!(prob_linear, solver)
    # solve!(prob_unreg, solver)
    # # Using functions in module:
    var_linear = Variable(1)
    v_linear_l1 = InsarTimeseries.invert_pixel(unw_vals, Blin, var_linear)
    var_unreg = Variable(size(B, 2))
    v_unreg_l1 = InsarTimeseries.invert_pixel(unw_vals, B, var_unreg)

    # v_linear_l1 = InsarTimeseries.invert_pixel(unw_vals, Blin, rho=1.0, alpha=1.5)
    # v_unreg_l1 = InsarTimeseries.invert_pixel(unw_vals, Blin, rho=1.0, alpha=1.5)

    return v_linear_lstsq, v_unreg_lstsq, v_linear_l1, v_unreg_l1
end

function integrate_velos(v_linear_lstsq, v_unreg_lstsq, v_linear_l1, v_unreg_l1)
    linear_lstsq = PHASE_TO_CM .* InsarTimeseries.integrate_velocities(v_linear_lstsq, timediffs)
    unreg_lstsq = PHASE_TO_CM .* InsarTimeseries.integrate_velocities(v_unreg_lstsq, timediffs)
    linear_l1 = PHASE_TO_CM .* InsarTimeseries.integrate_velocities(v_linear_l1, timediffs)
    unreg_l1 = PHASE_TO_CM .* InsarTimeseries.integrate_velocities(v_unreg_l1, timediffs)
    return linear_lstsq, unreg_lstsq, linear_l1, unreg_l1
end

function plot_insar(geolist, insar_linear_ts, insar_unreg_ts; title="", ylims=YLIMS)
    p = plot(geolist, insar_linear_ts, title=title, label="linear insar", ylabel="cm", 
             linewidth=3, ylims=ylims, legend=:bottomleft)
    plot!(p, geolist, insar_unreg_ts, label="unreg insar", marker=:o, linealpha=0.0)
    return p
end

function plot_gps!(p, dts, gps_los_data, gps_poly)
    day_nums = _get_day_nums(dts)
    plot!(p, dts, gps_los_data, marker=:x, color=:green, linealpha=0.0, markeralpha=0.5, markersize=1, label="gps")
    plot!(p, dts, gps_poly(day_nums), color=:green, linewidth=3, label="gps line fit")
end


function process_pixel(; station_name=nothing, plotting=false, reference_station=nothing)
    println("Processing $station_name")
    total_interval_days = (GEOLIST[end] - GEOLIST[1]).value

    # First solve for the velocities in L1 vs L2 to compare to GPS
    v_linear_lstsq, v_unreg_lstsq, v_linear_l1, v_unreg_l1 = solve_insar_ts(station_name)

    # dts, gps_los_data = get_gps_los(station_name, reference_station=REFERENCE_STATION)
    dts, gps_los_data = get_gps_los(station_name)
    # Convert to "days since start" for line fitting
    gps_poly = fit_line(dts, gps_los_data)
    offset, slope = Polynomials.coeffs(gps_poly)
    slope_gps_mm_yr = 365 * 10 * slope

    slope_insar_l2_mm_yr = (365 * 10 * PHASE_TO_CM * v_linear_lstsq)[1] # solution currently in phase
    slope_insar_l1_mm_yr = (365 * 10 * PHASE_TO_CM * v_linear_l1)[1]

    l2_error = slope_gps_mm_yr - slope_insar_l2_mm_yr
    l1_error = slope_gps_mm_yr - slope_insar_l1_mm_yr


    println("GPS velocity (mm / year): $slope_gps_mm_yr")
    println("InSAR L2 linear velocity (mm / year): $slope_insar_l2_mm_yr")
    println("InSAR L1 linear velocity (mm / year): $slope_insar_l1_mm_yr")
    println("==="^10)
    println("Difference (mm / yr): GPS - Insar L2 solution: $l2_error")
    println("Difference: (mm / yr) GPS - Insar L1 solution: $l1_error")
    println("==="^10)

    # Now get the full time series of insar by integrating
    linear_lstsq, unreg_lstsq, linear_l1, unreg_l1 = integrate_velos(v_linear_lstsq,
                                                                     v_unreg_lstsq,
                                                                     v_linear_l1,
                                                                     v_unreg_l1)

    if plotting
        # Plot the solutions vs gps
        p1 = plot_insar(GEOLIST, linear_lstsq, unreg_lstsq, title="L2 least squares solution")
        plot_gps!(p1, dts, gps_los_data, gps_poly)

        p2 = plot_insar(GEOLIST, linear_l1, unreg_l1, title="L1 norm minimization")
        plot_gps!(p2, dts, gps_los_data, gps_poly)

        plot(p1, p2)
        println("Saving $station_name.png")
        png(station_name)
    end
    return l1_error, l2_error
end


l1_errors, l2_errors = [], []

plotting = false

# station_name = station_name_list[1]
# station_name = "TXOZ"
# station_name_list = ["TXOZ"]

# l1_error, l2_error = process_pixel(station_name, UNW_STACK_FILE, GEOLIST, INTLIST, VALID_IGRAM_INDICES, plotting=plotting)
# append!(l2_errors, l2_error)
# append!(l1_errors, l1_error)
#
for station_name in station_name_list
    l1_error, l2_error = process_pixel(station_name=station_name, plotting=plotting)
    append!(l2_errors, l2_error)
    append!(l1_errors, l1_error)
end
println("TOTAL ERRORS (all in mm / year of velocity):")
println("L2 RMS error for all stations: $(rms(l2_errors))")
println("L1 RMS error for all stations: $(rms(l1_errors))")
println("L2 sum of abs errors for all stations: $(total_abs_error(l2_errors))")
println("L1 sum of abs errors for all stations: $(total_abs_error(l1_errors))")
println("L2 maximum errors : $(maximum(l2_errors))")
println("L1 maximum errors : $(maximum(l1_errors))")



# function check_convex_l2(Blin, unw_vals)
#     # To check that l2 norm min matches the backslash result
#     v2linear_test = Variable(1)
#     prob_linear_test2 = minimize(norm(Blin*v2linear_test - unw_vals, 2))
#     solve!(prob_linear_test2, ECOSSolver())
# end

# # The values already solved for:
# ts_linear = h5read(DEFO_FILENAME_LINEAR, "stack", (col, row, :))[1, 1, :];
# ts_unreg = h5read(DEFO_FILENAME, "stack", (col, row, :))[1, 1, :];
# p1 = plot(GEOLIST, ts_linear)
# plot!(p1, GEOLIST, ts_unreg)

