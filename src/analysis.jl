unw_vals_by_date(date, intlist, unw_vals) = unw_vals[date in intlist]

Blins_by_date(date, intlist, Blin) = Blin[date in intlist, :]

function remove_dates(bad_dates::Union{Date, Array{Date}}, intlist, unw_vals, B)
    good_idxs = _good_idxs(bad_dates, intlist)
    B_clean = B[good_idxs, :]
    unw_clean = unw_vals[good_idxs]
    return B_clean, unw_clean
end

_good_idxs(bad_date::Date, intlist) = .!(bad_date in intlist)

function _good_idxs(date_arr::Array{Date}, intlist)
    return .!reduce((x, y)-> x .| y, 
                    (b in intlist for b in date_arr),
                    init=falses(size(intlist)));
end

function solve_without_date(bad_dates::Union{Date, Array{Date}}, intlist, unw_vals, B)
    B_clean, unw_clean = remove_dates(bad_dates, intlist, unw_vals, B)
    velo_l1 = InsarTimeseries.invert_pixel(unw_clean, B_clean, rho=1.0, alpha=1.5)
    velo_lstsq = B_clean \ unw_clean
    # Return soluyion in mm/year
    return PHASE_TO_CM * 10 * 365 .* (velo_l1[1], velo_lstsq[1])
end

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

mean_abs_val(geolist, intlist, unw_vals) = [mean(abs.(unw_vals_by_date(d, intlist, unw_vals)))
                                             for d in geolist];
max_abs_val(geolist, intlist, unw_vals) = [maximum(abs.(unw_vals_by_date(d, intlist, unw_vals)))
                                             for d in geolist];


# Useful:
function remove_with_cutoff(diff_cutoff, geolist, intlist, unw_vals, B)
    l1_diffs, lstsq_diffs = compare_solutions(geolist, intlist, unw_vals, B)
    big_l1_diffs = l1_diffs .> diff_cutoff
    bad_days = geolist[big_l1_diffs]
    return solve_without_date(bad_days, intlist, unw_vals, B)
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
