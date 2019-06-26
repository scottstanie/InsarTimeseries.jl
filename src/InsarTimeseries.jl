"""Module InsarTimeseries

Testing docstring
"""
module InsarTimeseries

using Dates

"""Finds the number of days between successive .geo files"""
function day_diffs(geolist::Array{Date, 1})
	[difference.value for difference in diff(geolist)]
end

Array

"""
	matrixA(geolist::Array{Date, 1}, intlist::Array{Tuple{Date, Date}, 1}) -> Array{Int, 2}

Takes the list of igram dates and builds the SBAS A matrix

Returns the incident-like matrix from the SBAS paper: A*phi = dphi
	Each row corresponds to an igram, each column to a .geo
	value will be -1 on the early (slave) igrams, +1 on later (master)
"""
function matrixA(geolist, intlist)
    # We take the first .geo to be time 0, leave out of matrix
    # Match on date (not time) to find indices
	geolist = geolist[2:end]

	M = length(intlist)  # Number of igrams, 
	N = length(geolist)
	A = zeros(Int8, M, N)
	for j in 1:M
		early_date, late_date = intlist[j]
		idx_early = findfirst(isequal(early_date), geolist)
		if ~isnothing(idx_early)  # The first SLC will not be in the matrix
			A[j, idx_early] = -1
		end
		idx_late = findfirst(isequal(late_date), geolist)
        A[j, idx_late] = 1

	end
	A
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


end # module


# FUTURE TEST:
# geolist = collect(Date(2019,1,1):Day(2):Date(2019,1,7))
# intlist = []
# 
# for i in 1:length(geolist)-1
# 	for j in i+1:length(geolist)
# 		push!(intlist,(geolist[i], geolist[j]))
# 	end
# end
#
# julia> geolist
# #4-element Array{Date,1}:
#  2019-01-01
#  2019-01-03
#  2019-01-05
#  2019-01-07
# julia> intlist
# 6-element Array{Any,1}:
#  (2019-01-01, 2019-01-03)
#  (2019-01-01, 2019-01-05)
#  (2019-01-01, 2019-01-07)
#  (2019-01-03, 2019-01-05)
# (2019-01-03, 2019-01-07)
# (2019-01-05, 2019-01-07)
#
# A = InsarTimeseries.matrixA(geolist, intlist)
# julia> A = InsarTimeseries.matrixA(geolist, intlist)
# 6×3 Array{Int8,2}:
#   1   0  0
#   0   1  0
#   0   0  1
#  -1   1  0
#  -1   0  1
#   0  -1  1
