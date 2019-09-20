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

"""Get the mean values of igrams, per date, always taking the igram to be  (date, other)
(i.e. the sign of the phase is reversed if igram is (other, date)

The allows you to see if a date has some consistent delay added that is too
large to give useful phase difference values

Also reveals trends: For example, the final dates at the end of linear subsidence 
will have means with higher phase (positive phase = ground subsidence), 
since most igrams will look like (sunk ground -> neutral ground)
"""
function oneway_val(geolist, intlist, unw_vals, func=median)
    out = Array{Float64, 1}()
    for g in geolist
        # Make a 1 if g is the first date, -1 if second, and 0 otherwise
        mults = [g == ifg[1] ? 1 : (g == ifg[2] ? -1 : 0) for ifg in intlist]
        vals = filter(!iszero, mults .* unw_vals)
        push!(out, func(vals))
    end
    return out
end

"""Simpler means by day: just abs. value of phases, either mean or median"""
mean_abs_val(geolist, intlist, unw_vals) = [mean(abs.(vals_by_date(d, intlist, unw_vals)))
                                            for d in geolist];
median_abs_val(geolist, intlist, unw_vals) = [median(abs.(vals_by_date(d, intlist, unw_vals)))
                                              for d in geolist];


"""Used for robust variance est. (as a cutoff for outliers)
See https://en.wikipedia.org/wiki/Robust_measures_of_scale#IQR_and_MAD"""
iqr(arr) = quantile(arr, .75) - quantile(arr, .25)

""" 3* the sigma valud as calculated using MAD"""
mednsigma(arr, n=3) = n * mad(arr, normalize=true)  

two_way_cutoff(arr, nsigma) = (min(0, median(arr) - mednsigma(arr, nsigma)),  # Make sure it's negative
                               median(arr) + mednsigma(arr, nsigma))

# TODO: which of these do I really need?
function two_way_inliers(arr, nsigma)
    low, high = two_way_cutoff(arr, nsigma)
    return low .< arr .< high
end
function two_way_outliers(arr, nsigma)
    low, high = two_way_cutoff(arr, nsigma)
    return (arr .< low) .| (arr .> high)
end


"""Remove all igrams corresponding to the date with the highest mean"""
function peel_largest_n(geo, int, val, B; n=1)
    dates_to_remove = largest_n_dates(geo, int, val, n)
    return _remove_largest(geo, int, val, B, dates)
end

"""Remove all igrams corresponding to the dates with means falling outside an `nsigma` interval"""
function peel_nsigma(geo, int, val, B; nsigma=3)
    dates_to_remove = nsigma_days(geo, int, val, nsigma)
    return _remove_largest(geo, int, val, B, dates_to_remove)
end

"""Return the days of `geo` which are more than nsigma away from mean"""
function nsigma_days(geo, int, val, nsigma=3)
    # means = mean_abs_val(geo, int, val)
    # means = median_abs_val(geo, int, val); println("median abs")
    # @show mednsigma(means, nsigma)
    # TODO: change "means" if its not means
    # means = abs.(oneway_val(geo, int, val, mean))
    # means = abs.(oneway_val(geo, int, val, median))
    # means = oneway_val(geo, int, val, mean)
    
    means = oneway_val(geo, int, val, median); # println("median oneway")
    # @show mednsigma(means, nsigma)
    # TODO: change "means" if its not means

    # # FOR PRINTING ONLY
    # low, high = two_way_cutoff(means, nsigma)
    # println("Using cutoff around $(median(means)), spread $(mednsigma(means, nsigma)) : ($low, $high) ")

    return geo[two_way_outliers(means, nsigma)]
end

function _remove_largest(geo, int, val, B, dates_to_remove)
    int2, val2, B2 = remove_igrams(int, val, B, dates_to_remove)
    geo2 = [g for g in geo if !(g in dates_to_remove)]
    return geo2, int2, val2, B2
end


function prune_igrams(geolist, intlist, unw_pixel, B;
                      cor_pixel=nothing, cor_thresh=0.0,
                      mean_sigma_cutoff=3, 
                      fast_cm_cutoff=SENTINEL_WAVELENGTH/4)
    # 3 Reasons for removing igrams:
    # 1. remove all from .geo dates with huge averages (whole day is garbage)
    # 2. remove very low correlation
    #   NOTE: for now, we are skipping... doesn't seem to make much difference generally
    # 3. remove longer igrams w/ longer baseline than is possible
    #   e.g. if we have 2.5 cm/year velocity, anything longer than ~1 year
    #   will be all noise- we can't reliably sense such quickly moving ground
    # TODO: verify this third criteria?

    # @show "start", size(intlist)
    # 1. find outliers in this pixels' values and remove them
    geo_clean, intlist_clean, unw_clean, B_clean  = peel_nsigma(geolist, intlist, unw_pixel, B, nsigma=mean_sigma_cutoff)
    # @show "outlier: ", size(intlist_clean)

    # 2. low correlation cleaning
    if cor_thresh > 0 && !isnothing(cor_pixel)
        low_cor_igrams = intlist[cor_pixel .< cor_thresh]
        intlist_clean, unw_clean, B_clean  = remove_igrams(intlist_clean, unw_clean, B_clean, low_cor_igrams)
    end
    # @show "cor: ", size(intlist_clean)

    # 3. with rought velocity estimate, find igrams with too long of baseline
    # Here we assume beyond some wavelength fraction is too long to sense in one igram
    if fast_cm_cutoff < 5
        velo_cm = PHASE_TO_CM * (B_clean \ unw_clean)[1]  # cm / day
        # rho, alpha, abstol = 1.0, 1.6, 1e-3
        # velo_cm = PHASE_TO_CM * invert_pixel(unw_clean, B_clean, rho=rho, alpha=alpha, abstol=abstol)[1]
        day_cutoff = fast_cm_cutoff / abs(velo_cm)
        # @show (365*velo_cm), fast_cm_cutoff, day_cutoff
        too_long_igrams = [ig for ig in intlist_clean 
                           if temporal_baseline(ig) > day_cutoff]

        intlist_clean, unw_clean, B_clean  = remove_igrams(intlist_clean, unw_clean, B_clean, too_long_igrams)
    end
    # @show "fast: ", size(intlist_clean)
    return geo_clean, intlist_clean, unw_clean, B_clean
end


# OLDER WAY: (redundant code? TODO: cleanup)
"""Look for outliers in how much the solution shifts by just having a large mean value
Returns a Bool array the size of `geolist` with `true` marking the outliers"""
function find_mean_outliers(geolist, intlist, unw_vals, nsigma=3; B=nothing)
    # means = mean_abs_val(geolist, intlist, unw_vals)
    # means = abs.(oneway_val(geolist, intlist, unw_vals, mean))
    # means = mean_oneway_val(geolist, intlist, unw_vals)
    means = oneway_val(geolist, intlist, unw_vals, median)

    cutoff_val = median(means) + mednsigma(means, nsigma)
    println("Using cutoff of $cutoff_val: $(median(means)) + $(mednsigma(means, nsigma))")
    bad_idxs = means .> cutoff_val
    if iters > 1
        g2 = geolist[.!bad_idxs]
        i2, u2, B2 = remove_igrams(intlist, unw_vals, B, geolist[bad_idxs])
        # TODO: fix this
        return find_mean_outliers(g2, i2, u2, B=B2, iters=iters-1)
    end
    return bad_idxs
end

# function largest_n_dates(geo, int, val, n=length(geo))
#     # means = mean_abs_val(geo, int, val)
#     means = abs.(oneway_val(geo, int, unw_vals, mean))
# 
#     # Sort by the means to order the geolist
#     sorted_top = sort(collect(zip(means, geo)), rev=true)[1:n]
#     # Now upzip to get just the dates (and discard their means
#     return collect(zip(sorted_top...))[2]
# end

function solve_after_cutoff(geolist, intlist, unw_vals, B, nsigma=3; in_mm_yr=true)  #, direction=:high, method=:mean)
    bad_days = find_mean_outliers(geolist, intlist, unw_vals, nsigma)


    bad_days = geolist[bad_idxs]
    println("Removing $(length(bad_days)) days out of $(length(geolist)): $bad_days")
    return solve_without(bad_days, intlist, unw_vals, B, in_mm_yr=in_mm_yr)
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
# Outlier Removal through fancy means (loF = sklearn.neighbors.LOF
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

