import Polynomials
using Statistics: quantile
using StatsBase: mad

"""Get the values (unw, cc, etc.) of one date"""
vals_by_date(date::Date, intlist, vals) = vals[date in intlist]

"""Get the values (unw, cc, etc.) grouped by each date"""
vals_by_date(date_arr::Array{Date}, intlist, vals) = [vals[d in intlist] for d in date_arr]

Blins_by_date(date, intlist, Blin) = Blin[date in intlist, :]

function remove_igrams(intlist, unw_vals, B,
                       bad_items::Union{Date, AbstractArray{Date}, AbstractArray{Igram}}...)
    good_idxs = trues(size(intlist))
    for item in bad_items
        good_idxs .&= _good_idxs(item, intlist)
    end
    return remove_by_idx(good_idxs, intlist, unw_vals, B)
end


"""Pass directly the indexes of good ones igrams (removes all BUT the good_idxs)"""
function remove_by_idx(good_idxs::AbstractArray, intlist, unw_vals, B)
    B_clean = B[good_idxs, :]
    unw_clean = unw_vals[good_idxs]
    intlist_clean = intlist[good_idxs]
    return intlist_clean, unw_clean, B_clean 
end

_good_idxs(bad_date::Date, intlist) = .!(bad_date in intlist)

function _good_idxs(bad_date_arr::AbstractArray{Date}, intlist)
    return .!reduce((x, y)-> x .| y, 
                    (b in intlist for b in bad_date_arr),
                    init=falses(size(intlist)));
end
function _good_idxs(bad_igram_arr::AbstractArray{Igram}, intlist)
    return [!(ig in bad_igram_arr) for ig in intlist]
end

"""Solve without a geo date, array of dates, or array of igrams"""
function solve_without(bad_items::Union{Date, Array{Date}, Array{Igram}}, intlist, unw_vals, B; in_mm_yr=true)
    intlist_clean, unw_clean, B_clean  = remove_igrams(intlist, unw_vals, B, bad_items)
    velo_l1 = InsarTimeseries.invert_pixel(unw_clean, B_clean, rho=1.0, alpha=1.5)
    velo_lstsq = B_clean \ unw_clean
    # Return soluyion in mm/year if specified
    scale = in_mm_yr ? P2MM : 1

    return scale .* (velo_l1[1], velo_lstsq[1])
end


"""For one set of unw_vals, loop through each day of geolist and remove.
Solves once with no removals, and returns the difference of no removing
and removing each day. 
"""
function compare_solutions(geolist, intlist, unw_vals, B)
    base_l1, base_lstsq = solve_without(Date(2000,1,1), intlist, unw_vals, B)
    l1_diffs = Array{Float32, 1}(undef, length(geolist))
    lstsq_diffs = similar(l1_diffs)
    for (idx, d) in enumerate(geolist)
        l1, lstsq = solve_without(d, intlist, unw_vals, B)
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
    base_l1, base_lstsq = solve_without(Date(2000,1,1), intlist, unw_vals, B)
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

mean_abs_val(geolist, intlist, unw_vals) = [mean(abs.(vals_by_date(d, intlist, unw_vals)))
                                             for d in geolist];
max_abs_val(geolist, intlist, unw_vals) = [maximum(abs.(vals_by_date(d, intlist, unw_vals)))
                                             for d in geolist];

"""Used for robust variance est. (as a cutoff for outliers)
See https://en.wikipedia.org/wiki/Robust_measures_of_scale#IQR_and_MAD"""
iqr(arr) = quantile(arr, .75) - quantile(arr, .25)

""" 3* the sigma valud as calculated using MAD"""
mednsigma(arr, n=3) = n * mad(arr, normalize=true)  

"""Look for outliers in how much the solution shifts by just having a large mean value
Returns a Bool array the size of `geolist` with `true` marking the outliers"""
function find_mean_outliers(geolist, intlist, unw_vals, cutoff=3)
    means = mean_abs_val(geolist, intlist, unw_vals)

    cutoff_val = median(means) + mednsigma(means, cutoff)
    # println("Using cutoff of $cutoff")
    return means .> cutoff_val
end

function solve_after_cutoff(geolist, intlist, unw_vals, B, cutoff=3; in_mm_yr=true)  #, direction=:high, method=:mean)
    # if method == :mean
    bad_idxs = find_mean_outliers(geolist, intlist, unw_vals, cutoff)
    # elseif method == :diff
    #     bad_idxs = _remove_by_diff(geolist, intlist, unw_vals, B, cutoff)
    # end

    bad_days = geolist[bad_idxs]
    println("Removing $(length(bad_days)) days out of $(length(bad_idxs)): $bad_days")
    return solve_without(bad_days, intlist, unw_vals, B, in_mm_yr=in_mm_yr)
end

"""Remove all igrams corresponding to the date with the highest mean"""
function peel_largest_dates(geo, int, val, B, n=1)
    means = mean_abs_val(geo, int, val)
    dates_to_remove = largest_n_dates(geo, int, val, n)
    int2, val2, B2 = remove_igrams(int, val, B, dates_to_remove)
    geo2 = [g for g in geo if g != dates_to_remove]
    return geo2, int2, val2, B2
end

function largest_n_dates(geo, int, val, n=length(geo))
    means = mean_abs_val(geo, int, val)
    # Sort by the means to order the geolist
    sorted_top = sort(collect(zip(means, geolist)), rev=true)[1:n]
    # Now upzip to get just the dates (and discard their means
    return collect(zip(sorted_top...))[2]
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


######################
# Outlier Removal
######################
function find_outliers(B, vals, lof, n_neighbors=100, is1d=true)
    # lof = LOF(n_neighbors)
    X_train = is1d ? reshape(vals, :, 1) : cat(B, vals, dims=2)
    isout = _calc_outliers(X_train, lof)
    return B[isout, :], vals[isout]
end
_calc_outliers(X, lof) = lof.fit_predict(X) .== -1

function find_inliers(B, vals, lof, n_neighbors=100, is1d=true)
    # lof = LOF(n_neighbors)
    X_train = is1d ? reshape(vals, :, 1) : cat(B, vals, dims=2)
    isout = _calc_outliers(X_train, lof)
    return B[.!isout, :], vals[.!isout]
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

