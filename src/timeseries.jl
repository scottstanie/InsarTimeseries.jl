using Dates
using Statistics

const SENTINEL_WAVELENGTH = 5.5465763  # cm
const PHASE_TO_CM = SENTINEL_WAVELENGTH / (-4 * π )

"""Finds the number of days between successive .geo files"""
function day_diffs(geolist::Array{Date, 1})
	[difference.value for difference in diff(geolist)]
end



"""
	build_A_matrix(geolist::Array{Date, 1}, intlist::Array{Tuple{Date, Date}, 1}) -> Array{Int, 2}

Takes the list of igram dates and builds the SBAS A matrix

Returns the incident-like matrix from the SBAS paper: A*phi = dphi
	Each row corresponds to an igram, each column to a .geo
	value will be -1 on the early (slave) igrams, +1 on later (master)
"""
function build_A_matrix(geolist::Array{Date, 1}, intlist::Array{Tuple{Date, Date}, 1})
    # We take the first .geo to be time 0, leave out of matrix
    # Match on date (not time) to find indices
	geolist = geolist[2:end]

	M = length(intlist)  # Number of igrams, 
	N = length(geolist)
	A = zeros(Int8, M, N)
	for i in 1:M
		early_date, late_date = intlist[i]
		idx_early = findfirst(isequal(early_date), geolist)
		if ~isnothing(idx_early)  # The first SLC will not be in the matrix
			A[i, idx_early] = -1
		end
		idx_late = findfirst(isequal(late_date), geolist)
        A[i, idx_late] = 1

	end
	return A
end

"""Takes the list of igram dates and builds the SBAS B (velocity coeff) matrix

Each row corresponds to an igram, each column to a .geo
Values will be t_k+1 - t_k for columns after the -1 in A,
up to and including the +1 entry
"""
function build_B_matrix(geolist::Array{Date, 1}, intlist::Array{Tuple{Date, Date}, 1})

    timediffs = day_diffs(geolist)

    A = build_A_matrix(geolist, intlist)
	B = zeros(size(A))

	for i in 1:size(B, 1)
		row = A[i, :]
		start_idx = findfirst(row .== -1)
        # if no -1 entry, start at index 1. Otherwise
		if isnothing(start_idx)
			start_idx = 1
		else
			start_idx += 1
		end
		# +1 will always exist in row
        # End index is inclusive of the +1
		end_idx = findfirst(row .== 1)

        # Now fill in the time diffs in the index range
        B[i, start_idx:end_idx] = timediffs[start_idx:end_idx]
	end

    return B
end

"""Takes SBAS velocity output and finds phases

Arguments:
	velocities come from invert_sbas
	timediffs are the days between each SAR acquisitions
		length will be 1 less than num SAR acquisitions
"""
function integrate_velocities(velocities::Array{Float64, 1}, timediffs::Array{Int, 1})
    # multiply each column of vel array: each col is a separate solution
    phi_diffs = velocities .* td_ints

    # Now the final phase results are the cumulative sum of delta phis
	# This is equivalent to replacing missing with 0s (like np.ma.cumsum does)
	phi_arr = cumsum(coalesce.(phi_diffs, 0))

    # Add 0 as first entry of phase array to match geolist length on each col
	pushfirst!(phi_arr, 0)

    phi_arr
end


"""Subtracts reference pixel group from each layer

window is the size around the reference pixel to average at each layer
"""
function shift_stack(stack::Array{Float32, 3}, ref_row::Int, ref_col::Int; window::Int=3)
	half_win = div(window, 2)
	winsize = (2half_win + 1, 2half_win + 1)  # Make sure it is odd sized
	# temp = Array{Float32, 2}(undef, winsize)
	for k = 1:size(stack, 3)
		layer = view(stack, :, :, k)
		
		lil = layer[ref_row - half_win:ref_row + half_win, 
								ref_col - half_win:ref_col + half_win]
		layer_mean = mean(lil)
		if k == 1
			@show lil, layer_mean
		end
		stack[:, :, k] .-= layer_mean

	end
    return stack
end




