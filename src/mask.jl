function pruneA(A::Array{<:Integer, 2}, rows_to_delete)
    row_masks = trues(size(A, 1))
    row_masks[rows_to_delete] .= false
    A_short = A[row_masks, :]
end

function find_zero_columns(A_short)
    # Now find any oclumns that are all zero to delete
    column_sums = reshape(sum(abs2, A_short, dims=1), :)
    # bad_cols = findall(iszero.(column_sums))
    zero_columns = iszero.(column_sums)
end

function collapse_timediffs(timediffs, zero_columns)
    out_tds = Array{Int}(undef, 0)
    cur_td = 0
    for (idx, td) in enumerate(timediffs)
        cur_td += td
        if !zero_columns[idx]
            append!(out_tds, cur_td)
            cur_td = 0
        end
    end
    return out_tds
end
