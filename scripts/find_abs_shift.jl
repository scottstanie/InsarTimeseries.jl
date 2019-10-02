push!(LOAD_PATH,joinpath(expanduser("~/repos/InsarTimeseries.jl/src/")))
# import InsarTimeseries
# using InsarTimeseries: PHASE_TO_CM, STACK_FLAT_SHIFTED_DSET, STACK_FLAT_DSET, STACK_DSET, CC_FILENAME, UNW_FILENAME
using HDF5
using PyCall
using Printf: @printf

gps = InsarTimeseries.gps
sario = InsarTimeseries.sario
latlon = pyimport("apertools.latlon")

# PATH 78 
station_name_list78 = ["NMHB", "TXAD", "TXBG", "TXBL", "TXCE", "TXFS", "TXKM", 
                       "TXL2", "TXMC", "TXMH", "TXOE", "TXOZ", "TXS3", "TXSO"]
#
# PATH 85
station_name_list85 = ["TXKM", "TXMH", "TXFS", "TXAL"]
# BAD station reasons: MDO1 (nothing 2014-2017), TXVH (started 2018), TXPC (ended 2017)

# Functions for error evaluation
rms(arr::AbstractArray{T, 1}) where {T <: Any} = sqrt(mean(arr.^2))
rms(arr::AbstractArray{T, 2}) where {T <: Any} = [rms(arr[i, :]) for i in 1:size(arr, 1)]
maxabs(x) = maximum(abs.(x))

function print_errors(errors, station_name_list)
    println("Name | error ")
    println("=============")
    for (err, name) in zip(errors, station_name_list)
        # println("$name RESULTS:")
        @printf("%s |  %.2f \n", name, err)
    end
    println("==============")
    println("RMS  | Max-abs")
    @printf("%.2f |  %.2f \n", rms(errors), maxabs(errors))
end

function _get_station_rowcol(station_name)
    dem_rsc = sario.load("dem.rsc")
    lon, lat = gps.station_lonlat(station_name)
    return map(x-> convert(Int, x), 
               latlon.nearest_pixel(dem_rsc, lon=lon, lat=lat))
end


"""Load the insar data in a small patch around the pixel"""
function _get_val_at_station(insar_fname, station_name; dset="velos", window=5, kwargs...)
    # Note: swapping row and col due to julia/hdf5 transposes
    row, col = _get_station_rowcol(station_name)
    halfwin = div(window, 2)
    patch = h5read(insar_fname, dset, (col-halfwin:col+halfwin, row-halfwin:row+halfwin))
    return mean(patch)
end

get_station_values(fname, station_list::Array{String}; kwargs...) = [_get_val_at_station(fname, stat; kwargs...)
                                                                     for stat in station_list]

"""Function to take a velocity file and calculate the GPS errors at one station"""
function get_gps_error(insar_fname, station_name; dset="velos", window=5, ref_station=nothing, verbose=false)

    # insar derived solution in a small patch
    slope_insar_mm_yr = _get_val_at_station(insar_fname, station_name, dset=dset, window=window)

    slope_gps_mm_yr = InsarTimeseries.solve_gps_ts(station_name, ref_station)

    if verbose
        @show station_name
        println("GPS velocity (mm / year): $slope_gps_mm_yr")
        println("InSAR linear velocity (mm / year): $slope_insar_mm_yr")
        println("==="^10)
    end

    return slope_insar_mm_yr - slope_gps_mm_yr 
end
get_gps_error(fname, station_list::Array{String}; kwargs...) = [get_gps_error(fname, stat; kwargs...)
                                                                for stat in station_list]

"""Given a list of errors from insar-gps, find the const to add
to the insar solution image to minimize these errors
(converts the gps from relative to absolute)"""
function minimize_errors(error_list, search_range=-5:.1:5)
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

