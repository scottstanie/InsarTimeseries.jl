import Polynomials

unw_vals_by_date(date, intlist, unw_vals) = unw_vals[date in intlist]

Blins_by_date(date, intlist, Blin) = Blin[date in intlist, :]

function remove_dates(bad_dates::Union{Date, AbstractArray{Date}}, intlist, unw_vals, B)
    good_idxs = _good_idxs(bad_dates, intlist)
    return remove_dates(good_idxs, intlist, unw_vals, B)
end

"""Pass directly the indexes of good ones igrams (removes all BUT the good_idxs)"""
function remove_dates(good_idxs::AbstractArray, intlist, unw_vals, B)
    B_clean = B[good_idxs, :]
    unw_clean = unw_vals[good_idxs]
    intlist_clean = intlist_vals[good_idxs]
    return intlist_clean, unw_clean, B_clean 
end

_good_idxs(bad_date::Date, intlist) = .!(bad_date in intlist)

function _good_idxs(date_arr::Array{Date}, intlist)
    return .!reduce((x, y)-> x .| y, 
                    (b in intlist for b in date_arr),
                    init=falses(size(intlist)));
end

function solve_without_date(bad_dates::Union{Date, Array{Date}}, intlist, unw_vals, B)
    intlist_clean, unw_clean, B_clean  = remove_dates(bad_dates, intlist, unw_vals, B)
    velo_l1 = InsarTimeseries.invert_pixel(unw_clean, B_clean, rho=1.0, alpha=1.5)
    velo_lstsq = B_clean \ unw_clean
    # Return soluyion in mm/year
    return PHASE_TO_CM * 10 * 365 .* (velo_l1[1], velo_lstsq[1])
end


"""For one set of unw_vals, loop through each day of geolist and remove.
Solves once with no removals, and returns the difference of no removing
and removing each day. 
"""
function compare_solutions(geolist, intlist, unw_vals, B)
    base_l1, base_lstsq = solve_without_date(Date(2000,1,1), intlist, unw_vals, B)
    l1_diffs = Array{Float32, 1}(undef, length(geolist))
    lstsq_diffs = similar(l1_diffs)
    for (idx, d) in enumerate(geolist)
        l1, lstsq = solve_without_date(d, intlist, unw_vals, B)
        l1_diffs[idx] =  l1 - base_l1
        lstsq_diffs[idx] = lstsq - base_lstsq
    end
    return l1_diffs, lstsq_diffs
end

function compare_solutions_with_gps(geolist, intlist, unw_vals, station_name, linear=true)
    B = InsarTimeseries.build_B_matrix(geolist, intlist)
    B = linear ? sum(B, dims=2) : B

    l1_diffs = Array{Float32, 1}(undef, length(geolist))
    lstsq_diffs = similar(l1_diffs)

    slope_gps_mm_yr = solve_gps_ts(station_name, nothing)

    # First, solve with a dummy date so nothing is removed
    base_l1, base_lstsq = solve_without_date(Date(2000,1,1), intlist, unw_vals, B)
    base_l1_error = base_l1 - slope_gps_mm_yr
    base_lstsq_error = base_lstsq - slope_gps_mm_yr
    println("Base L1 error for $station_name = $base_l1_error")

    for (idx, d) in enumerate(geolist)
        l1, lstsq = solve_without_date(d, intlist, unw_vals, B)
        l1d = l1 - slope_gps_mm_yr
        lstsqd = lstsq - slope_gps_mm_yr

        # Note: if base is larger, the diff will be positive (an improvement)
        l1_diffs[idx] = abs(base_l1_error) - abs(l1d)
        lstsq_diffs[idx] = abs(base_lstsq_error) - abs(lstsqd)
    end
    return l1_diffs, lstsq_diffs, base_l1_error
end

mean_abs_val(geolist, intlist, unw_vals) = [mean(abs.(unw_vals_by_date(d, intlist, unw_vals)))
                                             for d in geolist];
max_abs_val(geolist, intlist, unw_vals) = [maximum(abs.(unw_vals_by_date(d, intlist, unw_vals)))
                                             for d in geolist];


# TODO Useful, but how to automaticall pick whether high,low, or abs??
function remove_with_cutoff(cutoff, geolist, intlist, unw_vals, B, direction=:high, method=:mean)
    if method == :mean
        bad_idxs = _remove_by_mean(geolist, intlist, unw_vals, cutoff)
    elseif method == :diff
        bad_idxs = _remove_by_diff(geolist, intlist, unw_vals, B, cutoff)
    end

    bad_days = geolist[bad_idxs]
    println("Removing $(length(bad_days)) days out of $(length(bad_idxs)): $bad_days")
    return solve_without_date(bad_days, intlist, unw_vals, B)
end

"""Look for outliers in how much the solution shifts by just having a large mean value"""
function _remove_by_mean(geolist, intlist, unw_vals, cutoff=nothing)
    means = mean_abs_val(geolist, intlist, unw_vals)
    cutoff = isnothing(cutoff) ? 1.5 * median(means) : cutoff
    return means .> cutoff
end

"""Look for outliers in how much the solution shifts by removing the single date"""
function _remove_by_diff(geolist, intlist, unw_vals, B, cutoff)
    l1_diffs, lstsq_diffs = compare_solutions(geolist, intlist, unw_vals, B)
    if direction == :high
        return lstsq_diffs .> cutoff
    elseif direction == :low
        return lstsq_diffs .< cutoff
    elseif direction == :abs
        return abs.(lstsq_diffs) .> cutoff
    else
        throw("direction must be :high, :low, or :abs")
    end
end


function geos_with_good_mean(geolist, intlist, unw_vals; cutoff=12)
    means_by_date = mean_abs_val(geolist, intlist, unw_vals)
    good_geo_idxs = [idx for (idx, val) in enumerate(means_by_date) if val < cutoff];
    bad_geo_idxs = [idx for (idx, val) in enumerate(means_by_date) if val >= cutoff];

    good_geos = geolist[good_geo_idxs]
    bad_geos = geolist[bad_geo_idxs]

    # good_igrams = filter(igram -> any(g in igram for g in good_geos), intlist)
    good_igrams = filter(igram -> !(igram in bad_geos), intlist)

    good_unw_idxs = indexin(good_igrams, intlist)
    good_unw_vals = unw_vals[good_unw_idxs]

    good_geos, good_igrams, good_unw_vals
end

import Base.-
-(x::Igram) = (x[2], x[1])

function phase_triplets(intlist, unw_vals)
    # triplets = Array{Tuple{Igram, eltype(unw_vals)}}(undef, ())
    triplets = Array{Tuple{Igram, Igram, Igram}, 1}()
    max_temp = maximum(temporal_baseline(intlist))
    for ii in 1:length(intlist)
        # val1 = unw_vals[ii]
        for jj in range(ii, stop=length(intlist))
            # val2 = unw_vals[jj]
            ig1, ig2 = intlist[ii], intlist[jj]
            third = (ig1[1], ig2[2])
            if (ig1[2] != ig2[1]) || temporal_baseline(third) > max_temp
                continue
            end
            # push!(triplets, third, vals)
            push!(triplets, (ig1, ig2, third))
        end
    end
    return triplets
end

############################
# GPS FUNCTIONS
############################
function get_gps_los(station_name, geo_path="../"; reference_station=nothing)
    dts, gps_los_data = InsarTimeseries.gps.load_gps_los_data(geo_path, station_name, 
                                                              start_year=2015, end_year=2018,
                                                              zero_mean=true, 
                                                              reference_station=reference_station)
    
    return [convert(Date, d) for d in dts], gps_los_data
end

"""Find the linear fit of MM per year of the gps station"""
function solve_gps_ts(station_name, reference_station=nothing)
    # NOTE: CURRENLT IGNORING THE REFERENCE STATION AND FORCING IT TO BE NOTHING
    dts, gps_los_data = get_gps_los(station_name, reference_station=nothing)
    # If we wanna compare with GPS subtracted too, do this:
    # dts, gps_los_data = get_gps_los(station_name, reference_station=reference_station)

    # Convert to "days since start" for line fitting
    gps_poly = fit_line(dts, gps_los_data)
    slope = length(gps_poly) == 2 ? Polynomials.coeffs(gps_poly)[2] : Polynomials.coeffs(gps_poly)[1]
    # offset, slope = Polynomials.coeffs(gps_poly)
    slope_gps_mm_yr = 365 * 10 * slope
end

_get_day_nums(dts) = [( d - dts[1]).value for d in dts]

function fit_line(dts, data)
    day_nums = _get_day_nums(dts)
    p = Polynomials.polyfit(day_nums, data, 1)
    # p(day_nums[end])
    return p
end

