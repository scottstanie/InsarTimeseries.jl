using InsarTimeseries
using Test



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
# A = InsarTimeseries.build_A_matrix(geolist, intlist)
# julia> A = InsarTimeseries.build_A_matrix(geolist, intlist)
# 6×3 Array{Int8,2}:
#   1   0  0
#   0   1  0
#   0   0  1
#  -1   1  0
#  -1   0  1
#   0  -1  1
