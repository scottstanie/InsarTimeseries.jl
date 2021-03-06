import Polynomials

Blins_by_date(date, intlist, Blin) = Blin[date in intlist, :]


"""Solve without a geo date, array of dates, or array of igrams"""
function solve_without(
    bad_items::Union{Date,Array{Date},Array{Igram}},
    geolist,
    intlist,
    unw_vals;
    in_mm_yr = true,
)
    intlist_clean, unw_clean = remove_igrams(intlist, unw_vals, bad_items)

    geolist_clean = [g for g in geolist if !(g in bad_items)]

    B_clean = sum(build_B_matrix(geolist_clean, intlist_clean), dims = 2)
    velo_l1 = invert_pixel_l1(unw_clean, B_clean)
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
    base_l1, base_lstsq = solve_without(Date(2000, 1, 1), geolist, intlist, unw_vals)
    l1_diffs = Array{Float32,1}(undef, length(geolist))
    lstsq_diffs = similar(l1_diffs)
    for (idx, d) in enumerate(geolist)
        l1, lstsq = solve_without(d, geolist, intlist, unw_vals)
        l1_diffs[idx] = l1 - base_l1
        lstsq_diffs[idx] = lstsq - base_lstsq
    end
    return l1_diffs, lstsq_diffs
end

"""Get the mean values of igrams, per date, always taking the igram to be  (date, other)
(i.e. the sign of the phase is reversed if igram is (other, date)

The allows you to see if a date has some consistent delay added that is too
large to give useful phase difference values

Also reveals trends: For example, the final dates at the end of linear subsidence 
will have means with higher phase (positive phase = ground subsidence), 
since most igrams will look like (sunk ground -> neutral ground)
"""
function oneway_val(geolist, intlist, unw_vals, func = median)
    out = Array{Float64,1}()
    for g in geolist
        # Make a 1 if g is the first date, -1 if second, and 0 otherwise
        mults = [g == ifg[1] ? 1 : (g == ifg[2] ? -1 : 0) for ifg in intlist]
        vals = filter(!iszero, mults .* unw_vals)
        push!(out, func(vals))
    end
    return out
end



# OLDER WAY: (redundant code? TODO: cleanup)
"""Look for outliers in how much the solution shifts by just having a large mean value
Returns a Bool array the size of `geolist` with `true` marking the outliers"""
function find_mean_outliers(geolist, intlist, unw_vals, nsigma = 3)
    # means = mean_abs_val(geolist, intlist, unw_vals)
    # means = abs.(oneway_val(geolist, intlist, unw_vals, mean))
    # means = mean_oneway_val(geolist, intlist, unw_vals)
    means = oneway_val(geolist, intlist, unw_vals, median)

    cutoff_val = median(means) + mednsigma(means, nsigma)
    println("Using cutoff of $cutoff_val: $(median(means)) + $(mednsigma(means, nsigma))")
    bad_idxs = means .> cutoff_val
    if iters > 1
        g2 = geolist[.!bad_idxs]
        i2, u2 = remove_igrams(intlist, unw_vals, geolist[bad_idxs])
        # TODO: fix this
        return find_mean_outliers(g2, i2, u2, iters = iters - 1)
    end
    return bad_idxs
end

"""Remove all igrams corresponding to the date with the highest mean"""
function peel_largest_n(geo, int, val; n = 1)
    dates_to_remove = largest_n_dates(geo, int, val, n)
    int2, val2 = remove_igrams(int, val, dates_to_remove)
    geo2 = [g for g in geo if !(g in dates_to_remove)]
    return geo2, int2, val2
end

function largest_n_dates(geo, int, val, n = length(geo))
    n < 1 && return Array{Date,1}()
    # means = mean_abs_val(geo, int, val)
    means = abs.(oneway_val(geo, int, val, median))
    println("largest abs(oneway)")

    # Sort by the means to order the geolist
    dates_tup = top_n(means, geo, n)[2]
    return [d for d in dates_tup]  # Array conversion
end

function top_n(values, keys, n, rev = true)
    sorted_top = sort(collect(zip(values, keys)), rev = rev)[1:n]
    # Now upzip to get just the dates (and discard their means
    return collect(zip(sorted_top...))
end

function solve_without_largest(geolist, intlist, unw_vals, B, n = 1; in_mm_yr = true)  # , direction=:high, method=:mean)
    bad_days = largest_n_dates(geolist, intlist, unw_vals, n)
    println("Removing $(length(bad_days)) days out of $(length(geolist)): $bad_days")
    return solve_without(bad_days, geolist, intlist, unw_vals, in_mm_yr = in_mm_yr)
end

function solve_after_cutoff(geolist, intlist, unw_vals, B, nsigma = 3; in_mm_yr = true)  # , direction=:high, method=:mean)
    bad_days = find_mean_outliers(geolist, intlist, unw_vals, nsigma)


    bad_days = geolist[bad_idxs]
    println("Removing $(length(bad_days)) days out of $(length(geolist)): $bad_days")
    return solve_without(bad_days, geolist, intlist, unw_vals, in_mm_yr = in_mm_yr)
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
    et = eltype(unw_vals)
    triplets = Array{Tuple{Igram,Igram,Igram},1}()
    # closures = Array{et,1}()
    phases = Array{Tuple{et,et,et},1}()
    # Note the longest one to avoid searching for ifgs that dont exist
    max_temp = maximum(temporal_baseline(intlist))
    for ii = 1:length(intlist)
        # val1 = unw_vals[ii]
        for jj in range(ii, stop = length(intlist))
            # val2 = unw_vals[jj]
            ig1, ig2 = intlist[ii], intlist[jj]
            third = (ig1[1], ig2[2])
            if (ig1[2] != ig2[1]) || temporal_baseline(third) > max_temp
                continue
            end
            kk = findfirst(isequal(third), intlist)
            val_tup = (unw_vals[ii], unw_vals[jj], -unw_vals[kk])
            any(v == zero(eltype(et)) for v in val_tup) && continue
            # push!(triplets, third, vals)
            push!(triplets, (ig1, ig2, third))
            push!(phases, val_tup)
            # push!(closures, unw_vals[ii] + unw_vals[jj] - unw_vals[kk])
        end
    end
    phase_arr = zeros(et, (length(triplets), 3))
    for (idx, val) in enumerate(phases)
        phase_arr[idx, :] .= val
    end
    closures = sum(phase_arr, dims = 2)
    return triplets, phase_arr, closures
end


######################
# Outlier Removal through fancy means (loF = sklearn.neighbors.LOF
######################
function find_outliers(B, vals, lof, n_neighbors = 100, is1d = true)
    # lof = LOF(n_neighbors)
    X_train = is1d ? reshape(vals, :, 1) : cat(B, vals, dims = 2)
    isout = _calc_outliers(X_train, lof)
    return B[isout, :], vals[isout]
end
_calc_outliers(X, lof) = lof.fit_predict(X) .== -1

function find_inliers(B, vals, lof, n_neighbors = 100, is1d = true)
    # lof = LOF(n_neighbors)
    X_train = is1d ? reshape(vals, :, 1) : cat(B, vals, dims = 2)
    isout = _calc_outliers(X_train, lof)
    return B[.!isout, :], vals[.!isout]
end
