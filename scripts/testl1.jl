push!(LOAD_PATH,joinpath(expanduser("~/repos/InsarTimeseries.jl/src/")))
# import InsarTimeseries
using InsarTimeseries: PHASE_TO_CM, STACK_FLAT_SHIFTED_DSET, STACK_FLAT_DSET, STACK_DSET
using Dates
using HDF5
using LinearAlgebra
using Statistics: mean
using Plots
using PyCall: pyimport
import Polynomials
import NaNMath
using Printf: @printf
pyplot()
# import PyPlot
# plt = PyPlot

nm = NaNMath
# import Convex
# import ECOS
# import SCS

p2mm = InsarTimeseries.PHASE_TO_CM * 365 * 10

# TODO: find the shift to match insar to gps best
# insar_arr_txmc = [365 * 10 * PHASE_TO_CM * solve_insar_ts(station_name, 5, "TXMC")[3][1] for station_name in  station_name_list]
# [ rms(gps_arr78 - (insar_arr_txmc .+ cc)) for cc in range(1, stop=3, length=10)]
#
YLIMS = (-12, 12)
# TODO: maybe get these names automatically?
# PATH 78 
station_name_list = ["NMHB", "TXAD", "TXBG", "TXBL", "TXCE", "TXFS", "TXKM", 
                     "TXL2", "TXMC", "TXMH", "TXOE", "TXOZ", "TXS3", "TXSO"]
#
# PATH 85
# station_name_list = ["TXKM", "TXMH", "TXFS", "TXAL"]
# BAD: MDO1 (nothing 2014-2017), TXVH (started 2018), TXPC (ended 2017)

# REFERENCE_STATION = "TXKM"
# REFERENCE_STATION = "TXAD"
REFERENCE_STATION = nothing

latlon = pyimport("apertools.latlon")
sario = pyimport("apertools.sario")
gps = pyimport("apertools.gps")

# const DEFO_FILENAME = "deformation_unreg_maxtemp400.h5"
# const DEFO_FILENAME_LINEAR = "deformation_linear_maxtemp400.h5"
# const DEFO_FILENAME_L1 = "deformation_linear_maxtemp400_l1.h5"
# println("Defo filename $DEFO_FILENAME, defo linear: $DEFO_FILENAME_LINEAR")

const UNW_STACK_FILE="unw_stack.h5"  # Or InsarTimeseries.UNW_FILENAME
const CC_STACK_FILE="cc_stack.h5"  # Or InsarTimeseries.CC_FILENAME
const ignore_geo_file = "geolist_ignore.txt"
max_temporal_baseline = 500
# max_temporal_baseline = 1200

GEOLIST, INTLIST, VALID_IGRAM_INDICES = InsarTimeseries.load_geolist_intlist(UNW_STACK_FILE, ignore_geo_file, max_temporal_baseline);
timediffs = InsarTimeseries.day_diffs(GEOLIST)


# For checking / plotting specific interestion points:
# Sinkhole
row, col = 1838, 1367
# row, col = [361, 257] .+ [1, 1]
# well uplift
row, col = 1223, 943
# row, col = [244, 188] .+ [1, 1]


function get_stack_vals(unw_stack_file::String, station_name::String, window=5,
                        dset=STACK_FLAT_DSET, valid_indices=VALID_IGRAM_INDICES;
                        reference_station=nothing)
    dem_rsc = sario.load("dem.rsc")
    lon, lat = gps.station_lonlat(station_name)
    row, col = map(x-> convert(Int, x), 
                   latlon.nearest_pixel(dem_rsc, lon=lon, lat=lat))

    # println("Getting data at $station_name: row $row, col $col, lon $lon, lat $lat")
    unw_vals = _read_vals(unw_stack_file, row, col, window, dset, valid_indices)
    return _subtract_reference(unw_vals, reference_station, unw_stack_file, window, dset, valid_indices)
end

function get_stack_vals(unw_stack_file::String, row::Int, col::Int, window=5,
                        dset=STACK_FLAT_DSET, valid_indices=VALID_IGRAM_INDICES;
                        reference_station=nothing)
    unw_vals = _read_vals(unw_stack_file, row, col, window, dset, valid_indices)
    return _subtract_reference(unw_vals, reference_station, unw_stack_file, window, dset, valid_indices)
end

function _read_vals(unw_stack_file, row, col, window=5, dset=STACK_FLAT_DSET, valid_indices=VALID_IGRAM_INDICES)
    # println("Loading $row, $col from $dset, avging window $window")
    halfwin = div(window, 2)
    # Note: swapping the col and row in h5read since julia is col-major
    unw_depth = h5read(unw_stack_file, dset, (col-halfwin:col+halfwin, row-halfwin:row+halfwin, :))
    unw_vals_all = vec(mean(unw_depth, dims=(1,2)))
    unw_vals = unw_vals_all[valid_indices]
end

function _subtract_reference(unw_vals, reference_station, unw_stack_file, window, dset, valid_indices)
    isnothing(reference_station) && return unw_vals
    return unw_vals - get_stack_vals(unw_stack_file, reference_station, window, dset, valid_indices)
end


function solve_insar_ts(station_name::String, window::Int=5, reference_station=nothing; cutoff=false)
    unw_vals = get_stack_vals(UNW_STACK_FILE, station_name, window, reference_station=reference_station)
    return solve_insar_ts(unw_vals, window, cutoff=cutoff)
end

function solve_insar_ts(row::Int, col::Int, window::Int=5, reference_station=nothing; cutoff=false)
    unw_vals = get_stack_vals(UNW_STACK_FILE, row, col, window, reference_station=reference_station)
    return solve_insar_ts(unw_vals, window, cutoff=cutoff)
end

function solve_insar_ts(unw_vals::Array{<:AbstractFloat, 1}, window::Int=5; cutoff=false)

    B = InsarTimeseries.build_B_matrix(GEOLIST, INTLIST)
    Blin = sum(B, dims=2)

    if cutoff
         v_linear_l1, v_linear_lstsq = InsarTimeseries.solve_after_cutoff(GEOLIST, INTLIST, unw_vals, Blin, in_mm_yr=false)
        return v_linear_lstsq, 0, v_linear_l1, 0
    end

    v_linear_lstsq = Blin \ unw_vals
    v_unreg_lstsq = B \ unw_vals

    # # solver = SCS.SCSSolver()
    # solver = ECOS.ECOSSolver(verbose=0)
    # # Using functions in module:
    # var_linear = Convex.Variable(1)
    # v_linear_l1 = InsarTimeseries.invert_pixel(unw_vals, Blin, var_linear)
    # var_unreg = Convex.Variable(size(B, 2))
    # v_unreg_l1 = InsarTimeseries.invert_pixel(unw_vals, B, var_unreg)

    # Using Huber Loss in ADMM
    v_linear_l1 = InsarTimeseries.invert_pixel(unw_vals, Blin, rho=1.0, alpha=1.5)
    v_unreg_l1 = InsarTimeseries.invert_pixel(unw_vals, Blin, rho=1.0, alpha=1.5)

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

_get_day_nums(dts) = [( d - dts[1]).value for d in dts]

function process_pixel(; station_name=nothing, plotting=false, reference_station=REFERENCE_STATION,
                       window=5, verbose=true, cutoff=false)
    total_interval_days = (GEOLIST[end] - GEOLIST[1]).value

    # First solve for the velocities in L1 vs L2 to compare to GPS
    reference_station = "TXKM"
    v_linear_lstsq, v_unreg_lstsq, v_linear_l1, v_unreg_l1 = solve_insar_ts(station_name, window, reference_station, cutoff=cutoff)

    # NOTE: CURRENLT IGNORING THE REFERENCE STATION AND FORCING IT TO BE NOTHING
    slope_gps_mm_yr = InsarTimeseries.solve_gps_ts(station_name, nothing)

    slope_insar_l2_mm_yr = (365 * 10 * PHASE_TO_CM * v_linear_lstsq)[1] # solution currently in phase
    slope_insar_l1_mm_yr = (365 * 10 * PHASE_TO_CM * v_linear_l1)[1]

    l2_error = slope_gps_mm_yr - slope_insar_l2_mm_yr
    l1_error = slope_gps_mm_yr - slope_insar_l1_mm_yr


    if verbose
        @show station_name reference_station
        println("GPS velocity (mm / year): $slope_gps_mm_yr")
        println("InSAR L2 linear velocity (mm / year): $slope_insar_l2_mm_yr")
        println("InSAR L1 linear velocity (mm / year): $slope_insar_l1_mm_yr")
        println("==="^10)
        println("Difference (mm / yr): GPS - Insar L2 solution: $l2_error")
        println("Difference: (mm / yr) GPS - Insar L1 solution: $l1_error")
        println("==="^10)
    end


    if plotting
        # Now get the full time series of insar by integrating
        linear_lstsq, unreg_lstsq, linear_l1, unreg_l1 = integrate_velos(v_linear_lstsq,
                                                                         v_unreg_lstsq,
                                                                         v_linear_l1,
                                                                         v_unreg_l1)
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

function calc_error_matrix(station_name_list, plotting=false)
    println("Calculating errors at all station/reference combos")
    l1_error_matrix = zeros(length(station_name_list), length(station_name_list))
    l2_error_matrix = zeros(length(station_name_list), length(station_name_list))

    for (i, ref_name) in enumerate(station_name_list)
        println("Calculating errors for $ref_name")
        for (j, station_name) in enumerate(station_name_list)
            # println("Processing $station_name")
            l1_error, l2_error = process_pixel(station_name=station_name, plotting=plotting, 
                                               reference_station=ref_name, verbose=false)
            l1_error_matrix[i, j] = l1_error
            l2_error_matrix[i, j] = l2_error
        end
    end
    return l1_error_matrix, l2_error_matrix
end

# Reminder: error matrix rows correspond to 1 single reference station
function print_matrix_stats(l1_error_matrix, l2_error_matrix, station_name_list)
    println("Errors compared to GPS when using <station name> as the reference point")
    for (mat, lname) in zip((l2_error_matrix, l1_error_matrix), ("L2", "L1"))
        println("$lname RESULTS:")
        println("Name | RMS error | Total abs | Max abs")
        for i in 1:length(station_name_list)
            row = mat[i, :]
            a, b, c, d = station_name_list[i], rms(row), total_abs_error(row), maximum(abs.(row))
            @printf("%s :    %.2f     |   %.2f   |    %.2f \n", a, b, c, d)
        end
    end
end

# Functions for error evaluation
rms(arr::AbstractArray{Any, 1}) = sqrt(mean(arr.^2))
rms(arr::AbstractArray{Any, 2}) = [rms(arr[i, :]) for i in 1:size(arr, 1)]
total_abs_error(arr::AbstractArray{Any, 1}) = sum(abs.(arr))
total_abs_error(arr::AbstractArray{Any, 2}) = sum(abs.(arr), dims=2)



plotting = false

############
# Try and remove a day, check all GPS errors
############

num_geos = length(GEOLIST)
num_stations = length(station_name_list)

l1_diff_matrix = zeros(num_stations, num_geos)
base_errors = zeros(num_stations)

for (idx, station_name) in enumerate(station_name_list)
    window = 5
    # I'm assuming TXKM is the best based on other analysis here
    ref_stat = "TXKM"
    # unw_vals = get_stack_vals(UNW_STACK_FILE, station_name, window, reference_station=ref_stat)

    # The "diff" being positive means an improvement (new error off GPS is lower), while
    # negative means the new err from GPS is bigger than before
    #
    # l1_diffs, l2_diffs, base_l1_error = InsarTimeseries.compare_solutions_with_gps(GEOLIST, INTLIST, unw_vals, station_name)

    # l1_diff_matrix[idx, :] = l1_diffs
    # base_errors[idx] = base_l1_error
end

# # top row is RMS, bottom row is max
# new_errs = zeros(2, num_geos)
# for idx in 1:num_geos
#     # if l1_diff_matrix entry is positive, it is an improvement, 
#     # so we subtract from base.
#     new_err = abs.(abs.(base_errors) .- l1_diff_matrix[:, idx])
#     # reminder: l1_diffs[idx] = abs(base_l1_error) - abs(l1d)
#     # which means we are getting back to `l1d = l1 - slope_gps_mm_yr`
# 
#     # println("$(rms(new_err)) , $(maximum(new_err))")
#     new_errs[:, idx] = [rms(new_err) , maximum(new_err)]
# end

# bases = [rms(abs.(base_errors)), maximum(abs.(base_errors))]
# err_diffs = new_errs .- bases;

function plot_errs(err_diffs, geolist)
    scatter(geolist, err_diffs[1, :], label="RMS", title="difference in errors (negative=improvement)")
    scatter!(geolist, err_diffs[2, :], label="maximum")
end


###### 

l1_errors, l2_errors = [], []

# Compute all combos:
# @time l1_error_matrix, l2_error_matrix = calc_error_matrix(station_name_list, plotting)
# print_matrix_stats(l1_error_matrix, l2_error_matrix, station_name_list)

# station_name = station_name_list[1]
# station_name = "TXOZ"
# station_name_list = ["TXMH"]

# l1_error, l2_error = process_pixel(station_name, UNW_STACK_FILE, GEOLIST, INTLIST, VALID_IGRAM_INDICES, plotting=plotting)
# append!(l2_errors, l2_error)
# append!(l1_errors, l1_error)

for station_name in station_name_list
    l1_error, l2_error = process_pixel(station_name=station_name, plotting=plotting, cutoff=true)
    append!(l2_errors, l2_error)
    append!(l1_errors, l1_error)
end

println("TOTAL ERRORS (all in mm / year of velocity):")
println("L2 RMS error for all stations: $(rms(l2_errors))")
println("L1 RMS error for all stations: $(rms(l1_errors))")
println("L2 sum of abs errors for all stations: $(total_abs_error(l2_errors))")
println("L1 sum of abs errors for all stations: $(total_abs_error(l1_errors))")
println("L2 maximum errors : $(maximum(abs.(l2_errors)))")
println("L1 maximum errors : $(maximum(abs.(l1_errors)))")



# function check_convex_l2(Blin, unw_vals)
#     # To check that l2 norm min matches the backslash result
#     v2linear_test = Convex.Variable(1)
#     prob_linear_test2 = minimize(norm(Blin*v2linear_test - unw_vals, 2))
#     Convex.solve!(prob_linear_test2, ECOSSolver())
# end

# # The values already solved for:
# ts_linear = h5read(DEFO_FILENAME_LINEAR, "stack", (col, row, :))[1, 1, :];
# ts_unreg = h5read(DEFO_FILENAME, "stack", (col, row, :))[1, 1, :];
# p1 = plot(GEOLIST, ts_linear)
# plot!(p1, GEOLIST, ts_unreg)

function plotunws(B1, B2, unw1, unw2)
    yh = maximum(vcat(unw1, unw2)) + 1
    yl = minimum(vcat(unw1, unw2)) - 1
    p1 = plot(B1, unw1, line=:nothing, marker=:o, color=:blue, ylims=(yl, yh))
    p2 = plot(B2, unw2, line=:nothing, marker=:o, color=:green, ylims=(yl, yh))
    plot(p1, p2)
end

function load_multi_temp(baselines...; station_name=nothing, rowcol=nothing, window=1)
    intlists, idxs, vals, Bs, ccvals = [], [], [], [], []
    for baseline in baselines
        geolist, cur_intlist, cur_valid_idxs = InsarTimeseries.load_geolist_intlist(UNW_STACK_FILE, ignore_geo_file, baseline);
        if isnothing(station_name)
            cur_unw_vals = get_stack_vals(UNW_STACK_FILE, rowcol..., window, STACK_FLAT_DSET,
                                        cur_valid_idxs, reference_station="TXKM");
            # TODO: add back in when the stack sizes match up again
            # ccval = get_stack_vals(UNW_STACK_FILE, rowcol..., window, STACK_DSET, cur_valid_idxs)
        else
            cur_unw_vals = get_stack_vals(UNW_STACK_FILE, station_name, window, STACK_FLAT_DSET,
                                        cur_valid_idxs, reference_station="TXKM");
            # TODO: add back in when the stack sizes match up again
            # ccval = get_stack_vals(UNW_STACK_FILE, station_name, window, STACK_DSET, cur_valid_idxs)
        end
        B = InsarTimeseries.build_B_matrix(geolist, cur_intlist)

        push!(intlists, cur_intlist)
        push!(idxs, cur_valid_idxs)
        push!(vals, cur_unw_vals)
        push!(Bs, sum(B, dims=2))
        # push!(ccvals, ccval)
    end
    return intlists, idxs, vals, Bs # , ccvals
end

function plot_multi_temp(Bs, vals, temps, title=".unw vals at different baselines")
    n = length(Bs)

    sort!(temps, rev=true)
    sort!(Bs, by=x->length(x), rev=true)
    sort!(vals, by=x->length(x), rev=true)

    ym = maximum(maximum(abs.(v)) for v in vals) + 3
    # fig, axes = plt.subplots(1, n)
    fig, axes = plt.subplots()
    for ii = 1:n
        # axes[ii].plot(Bs[ii], vals[ii], ".")
        axes.plot(Bs[ii], vals[ii], ".", label=temps[ii])
    end
    axes.set_ylim((-ym, ym))
    axes.legend()
    axes.set_title(title)
    return fig, axes
end
