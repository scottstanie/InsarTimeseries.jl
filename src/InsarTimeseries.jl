module InsarTimeseries

# import Pkg; Pkg.add("Dates")
using Dates


function matrixA(geolist, intlist)
    """Takes the list of igram dates and builds the SBAS A matrix

    Args:
	geolist (Array{Date}, 1): datetimes of the .geo acquisitions
	intlist (Array{Tuple{Date, Date}, 1}

    Returns:
		Array{Int, 2}
        np.array 2D: the incident-like matrix from the SBAS paper: A*phi = dphi
            Each row corresponds to an igram, each column to a .geo
            value will be -1 on the early (slave) igrams, +1 on later (master)
    """
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
end # matrixA



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
