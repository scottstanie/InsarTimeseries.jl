

# # TODO: is this nunecessary? Maybe i'll need once I start doing unregularized...
# """Takes the indices of igrams to ignore, then shrinks down the full size B matrix and
# timediffs array to remove zeros
# """
# function cut_igram_dates(A, B, timediffs, ignore_int_indices)
#     row_masks = find_row_masks(size(A, 1), ignore_int_indices)
#     column_masks = find_column_masks(A[row_masks, :])
# 
#     timediffs_cut = collapse_timediffs(timediffs, column_masks)
#     A_cut = A[row_masks, column_masks]
#     return timediffs_cut, B_cut
# end
# 
# function find_row_masks(num_rows, rows_to_delete)
#     row_masks = trues(num_rows)
#     row_masks[rows_to_delete] .= false
#     return row_masks
# end
# 
# 
# """Take the row-deleted A, find columns that lost all igrams to delete"""
# function find_column_masks(A_short)
#     # Now find any columns that are all zero to delete, make it 1D
#     column_sums = vec(sum(abs2, A_short, dims=1))
# 
#     # Trues are good columns (nonzero), false is a bad column
#     column_masks = .!iszero.(column_sums)
# end
# 
# 
# """Using the missing columns, recalculate the time diff vector"""
# function collapse_timediffs(timediffs, column_masks)
#     out_tds = Array{Int}(undef, 0)
#     cur_td = 0
#     for (idx, td) in enumerate(timediffs)
#         cur_td += td
#         if column_masks[idx]
#             append!(out_tds, cur_td)
#             cur_td = 0
#         end
#     end
#     return out_tds
# end
# 
