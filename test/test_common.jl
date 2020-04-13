using Dates

# Test data is 3 days after a starting reference on Jan 1:
# The actual phi phase is [0, 1, 2, 3]
#
#
geolist = [Date(2019, 1, 1), Date(2019, 1, 3), Date(2019, 1, 6), Date(2019, 1, 10)]
intlist = [(geolist[i], geolist[j]) for i = 1:length(geolist)-1 for j = i+1:length(geolist)]

# julia> intlist
#  (2019-01-01, 2019-01-03)  
#  (2019-01-01, 2019-01-06)
#  (2019-01-01, 2019-01-10)
#  (2019-01-03, 2019-01-06)
#  (2019-01-03, 2019-01-10)
#  (2019-01-06, 2019-01-10)

@test span(geolist) == 9
@test span(intlist) == 9
