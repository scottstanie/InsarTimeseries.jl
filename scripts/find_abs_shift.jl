push!(LOAD_PATH, joinpath(expanduser("~/repos/InsarTimeseries.jl/src/")))
# import InsarTimeseries
# using InsarTimeseries: PHASE_TO_CM, STACK_FLAT_SHIFTED_DSET, STACK_FLAT_DSET, STACK_DSET, CC_FILENAME, UNW_FILENAME
import Polynomials
import Glob
using HDF5
using PyCall
using Printf: @printf
try
    using MapImages
catch
end

# gps = InsarTimeseries.gps
latlon = pyimport("apertools.latlon")
gps = pyimport("apertools.gps")

reloadpython(modname) = pyimport("importlib")."reload"(modname)
reloadpython(gps)

LOS_MAP_FILE = "los_enu.tif"

# PATH 78 
station_name_list78 = [
    "NMHB",
    "TXAD",
    "TXBG",
    "TXBL",
    "TXCE",
    "TXFS",
    "TXKM",
    "TXL2",
    "TXMC",
    "TXMH",
    "TXOE",
    "TXOZ",
    "TXS3",
    "TXSO",
]
#
# PATH 85
station_name_list85 = ["TXKM", "TXMH", "TXFS", "NMHB", "TXAD", "TXS3"]
# Take out TXAL cuz it's not in crop
# BAD station reasons: MDO1 (nothing 2014-2017), TXVH (started 2018), TXPC (ended 2017)

all_stations = sort(unique([station_name_list78; station_name_list85]))

# Functions for error evaluation
rms(arr::AbstractArray{T,1}) where {T <: Any} = sqrt(mean(filter(!isnan, arr.^2)))
# rms(arr::AbstractArray{T, 2}) where {T <: Any} = [rms(arr[i, :]) for i in 1:size(arr, 1)]
maxabs(x) = maximum(abs.(filter(!isnan, x)))

# Function to take a velocity file and calculate the GPS errors at one station
function get_gps_error(
    insar_fname::String,
    station_name::String;
    dset="velos",
    window=5,
    ref_station=nothing,
    verbose=false,
    shift=0.0,
    negate=false,
)
    # If it's a file with LOS flipped:
    flip_sign = negate ? -1 : 1

    # insar derived solution in a small patch
    slope_insar_mm_yr =
        flip_sign * _get_val_at_station(insar_fname, station_name, dset=dset, window=window) + shift
    slope_insar_mm_yr -= isnothing(ref_station) ? 0 :
        (
        flip_sign * _get_val_at_station(insar_fname, ref_station, dset=dset, window=window) + shift
    )

    start_date, end_date = _get_date_range(insar_fname, dset)
    slope_gps_mm_yr = solve_gps_ts(
        station_name,
        ref_station,
        start_date=start_date,
        end_date=end_date,
        los_map_file=joinpath(dirname(insar_fname), LOS_MAP_FILE),
    )

    if verbose
        @show station_name
        println("GPS velocity (mm / year): $slope_gps_mm_yr")
        println("InSAR linear velocity (mm / year): $slope_insar_mm_yr")
        println("==="^10)
    end


    return slope_insar_mm_yr - slope_gps_mm_yr
end
get_gps_error(fname::String, station_list::Array{String}; kwargs...) =
    get_gps_error.(fname, station_list; kwargs...)
# ! note: dot broadcasting is same as [get_gps_error.(fname, stat; kwargs...) for stat in station_list]

_get_date_range(fname::AbstractString) = extrema(Sario.load_geolist_from_h5(fname))
_get_date_range(fname::AbstractString, dset::AbstractString) =
    extrema(Sario.load_geolist_from_h5(fname, dset))

function _get_station_rowcol(station_name, directory=".")
    demrsc = Sario.load(joinpath(directory, "dem.rsc"))
    return map(x -> convert(Int, x), MapImages.station_rowcol(station_name, demrsc))
end


"""Load the insar data in a small patch around the pixel"""
function _get_val_at_station(
    insar_fname,
    station_name;
    dset="velos",
    window=5,
    kwargs...,
)
    # Note: swapping row and col due to julia/hdf5 transposes
    row, col = _get_station_rowcol(station_name, dirname(insar_fname))
    halfwin = div(window, 2)
    try
        patch =
            h5read(insar_fname, dset, (col - halfwin:col + halfwin, row - halfwin:row + halfwin))
        return mean(patch)
    catch e
        println(e)
        return NaN
    end
end

get_station_values(fname, station_list::Array{String}; kwargs...) =
    [_get_val_at_station(fname, stat; kwargs...) for stat in station_list]

"""Given a list of errors from insar-gps, find the const to add
to the insar solution image to minimize these errors
(converts the gps from relative to absolute)"""
function minimize_errors(error_list, search_range=-5:0.1:5)
    best_rms, best_maxabs = 100, 100
    c_rms, c_maxabs = 0, 0
    for c in search_range
        cur_rms = rms(c .+ error_list)
        cur_maxabs = maxabs(c .+ error_list)
        if cur_rms < best_rms
            c_rms = c
            best_rms = cur_rms
        end

        if cur_maxabs < best_maxabs
            c_maxabs = c
            best_maxabs = cur_maxabs
        end
    end
    return c_rms, c_maxabs, best_rms, best_maxabs
end

function print_errors(
    errors78::AbstractArray{<:AbstractFloat},
    errors85::AbstractArray{<:AbstractFloat},
    station_name_list::AbstractArray{<:AbstractString},
)
    println("Name | error 78 | error 85 ")
    println("============================")
    for (err78, err85, name) in zip(errors78, errors85, station_name_list)
        # println("$name RESULTS:")
        @printf("%s |  %.2f | %.2f \n", name, err78, err85)
    end
    println("==============")
    println("RMS 78 | Max-abs 78 | RMS 85 | Max-abs 85")
    @printf(
        "%.2f |  %.2f | %.2f |  %.2f  \n",
        rms(errors78),
        maxabs(errors78),
        rms(errors85),
        maxabs(errors85)
    )
end


function compare_solutions_with_gps(geolist, intlist, unw_vals, station_name, linear=true)
    B = prepB(geolist, intlist, constant_velocity=linear)

    l1_diffs = Array{Float32,1}(undef, length(geolist))
    lstsq_diffs = similar(l1_diffs)

    slope_gps_mm_yr = solve_gps_ts(station_name, nothing)

    # First, solve with a dummy date so nothing is removed
    base_l1, base_lstsq = solve_without(Date(2000, 1, 1), intlist, unw_vals, B)
    base_l1_error = base_l1 - slope_gps_mm_yr
    base_lstsq_error = base_lstsq - slope_gps_mm_yr
    println("Base L1 error for $station_name = $base_l1_error")

    for (idx, d) in enumerate(geolist)
        l1, lstsq = solve_without(d, intlist, unw_vals, B)
        l1d = l1 - slope_gps_mm_yr
        lstsqd = lstsq - slope_gps_mm_yr

        # Note: if base is larger, the diff will be positive (an improvement)
        l1_diffs[idx] = abs(base_l1_error) - abs(l1d)
        lstsq_diffs[idx] = abs(base_lstsq_error) - abs(lstsqd)
    end
    return l1_diffs, lstsq_diffs, base_l1_error
end

############################
# GPS FUNCTIONS
############################
function get_gps_los(
    station_name;
    los_map_file=LOS_MAP_FILE,
    reference_station=nothing,
    start_date=Date(2014, 11, 1),
    end_date=Date(2019, 1, 1),
    days_smooth=0,
)
    dts, gps_los_data = gps.load_gps_los_data(
        days_smooth=days_smooth,
        los_map_file=los_map_file,
        station_name=station_name,
        start_date=start_date,
        end_date=end_date,
        zero_mean=true,
        reference_station=reference_station,
    )

    return [convert(Date, d) for d in dts], gps_los_data
end

function get_gps_enu(station_name)
    dts, enu_df = gps.load_station_enu(
        station_name,
        start_date=Date(2014, 11, 1),
        end_date=Date(2019, 1, 1),
        zero_mean=true,
    )

    # Convert from PyObjects to Arrays
    dts = [convert(Date, d) for d in dts]
    east = [r for r in enu_df.east]
    north = [r for r in enu_df.north]
    up = [r for r in enu_df.up]
    return dts, east, north, up
end


"""Find the linear fit of MM per year of the gps station"""
function solve_gps_ts(
    station_name,
    reference_station=nothing;
    start_date=Date(2014, 11, 1),
    end_date=Date(2019, 1, 1),
    los_map_file=LOS_MAP_FILE,
)
    dts, gps_los_data = get_gps_los(
        station_name,
        reference_station=reference_station,
        start_date=start_date,
        end_date=end_date,
        los_map_file=los_map_file,
    )
    # If we wanna compare with GPS subtracted too, do this:
    # dts, gps_los_data = get_gps_los(station_name, reference_station=reference_station)

    # Convert to "days since start" for line fitting
    gps_poly = fit_line(dts, gps_los_data)
    slope = length(gps_poly) == 2 ? Polynomials.coeffs(gps_poly)[2] :
        Polynomials.coeffs(gps_poly)[1]
    # offset, slope = Polynomials.coeffs(gps_poly)
    slope_gps_mm_yr = 365 * 10 * slope
    return slope_gps_mm_yr
end

function fit_line(dts, data)
    day_nums = InsarTimeseries._get_day_nums(dts)
    p = Polynomials.polyfit(day_nums, data, 1)
    # p(day_nums[end])
    return p
end

function shift_from(src::HDF5Dataset, shift::Real)
    tmp = src .+ shift
    tmp[src .== 0] .= 0
    return tmp
end

function shift_dset(shift, fname, src::String="velos/1", dest="velos_shifted/1")
    h5open(fname, "r+") do f
        out = shift_from(f[src], shift)
        f[dest] = out
    end
    h5writeattr(fname, dest, h5readattr(fname, src))
    return
end

function scale_dset!(scale, fname, src::String="velos/1")
    a = h5readattr(fname, src)
    h5open(fname, "cw") do f
        tmp = scale * read(f, src)
        o_delete(f, src)
        f[src] = tmp
    end
    h5writeattr(fname, src, a)
    return
end

function get_shift(fname, src::String="velos/1", dest="velos_shifted/1")
    h5open(fname, "r+") do f
        return f[dest][div(end, 2), div(end, 2)] - f[src][div(end, 2), div(end, 2)]
    end
    end

function save_as_unw(fname, dset="velos/1")
    amp = abs.(Sario.load(Glob.glob("*.int")[1]))
    img = Sario.load(fname, dset_name=dset)
    outname = replace(fname, ".h5" => ".unw")
    Sario.save(outname, cat(amp, Float32.(img ./ 10), dims=3))
end
