read_geolist_file(filename::String) = sario.find_geos(filename=filename)

"""Read extra file to ignore certain dates of interferograms"""
function find_valid_indices(geo_date_list::Array{Date, 1}, igram_date_list::Array{Igram, 1}, ignore_geo_file, max_temporal_baseline)
    if isnothing(ignore_geo_file) && isnothing(max_temporal_baseline)
        return 1:length(geo_date_list), 1:length(igram_date_list)
    end

    ignore_geos = sort(sario.find_geos(filename=ignore_geo_file, parse=true))
    println("Ignoring the following .geo dates:")
    println(ignore_geos)

    # First filter by remove igrams with either date in `ignore_geo_file`
    valid_geos = [g for g in geo_date_list if !(g in ignore_geos)]
    valid_igrams = [i for i in igram_date_list if !((i[1] in ignore_geos) || (i[2] in ignore_geos))]

    # Now also remove igrams spanning longer than `max_temporal_baseline`
    #
    # Compute for logging purposes:
    if !isnothing(max_temporal_baseline)
        too_long_igrams = filter(ig -> temporal_baseline(ig) > max_temporal_baseline, valid_igrams)
        println("Ignoring $(length(too_long_igrams)) igrams with longer baseline than $max_temporal_baseline days")
        valid_igrams = filter(ig -> temporal_baseline(ig) <= max_temporal_baseline, valid_igrams)
    end

    valid_geo_indices = indexin(valid_geos, geo_date_list)
    valid_igram_indices = indexin(valid_igrams, igram_date_list)

    println("Ignoring $(length(igram_date_list) - length(valid_igrams)) igrams total")

    return valid_geo_indices, valid_igram_indices
end


function find_row_masks(num_rows, rows_to_delete)
    row_masks = trues(num_rows)
    row_masks[rows_to_delete] .= false
    return row_masks
end


"""Take the row-deleted A, find columns that lost all igrams to delete"""
function find_column_masks(A_short)
    # Now find any columns that are all zero to delete
    column_sums = sum(abs2, A_short, dims=1)
    # Make sure it is 1D, not 2D
    column_sums = reshape(column_sums, :)

    # Trues are good columns (nonzero), false is a bad column
    column_masks = .!iszero.(column_sums)
end


"""Using the missing columns, recalculate the time diff vector"""
function collapse_timediffs(timediffs, column_masks)
    out_tds = Array{Int}(undef, 0)
    cur_td = 0
    for (idx, td) in enumerate(timediffs)
        cur_td += td
        if column_masks[idx]
            append!(out_tds, cur_td)
            cur_td = 0
        end
    end
    return out_tds
end

"""Takes the indices of igrams to ignore, then shrinks down the full size B matrix and
timediffs array to remove zeros
"""
function cut_igram_dates(A, B, timediffs, ignore_int_indices)
    row_masks = find_row_masks(size(A, 1), ignore_int_indices)
    column_masks = find_column_masks(A[row_masks, :])

    timediffs_cut = collapse_timediffs(timediffs, column_masks)
    A_cut = A[row_masks, column_masks]
    return timediffs_cut, B_cut
end

