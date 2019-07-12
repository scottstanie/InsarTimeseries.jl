using Dates

geolist = [Date(2019,1,1), Date(2019,1,3), Date(2019,1,6), Date(2019,1,10)]
intlist = [(geolist[i], geolist[j]) for i in 1:length(geolist)-1 for j in i+1:length(geolist)]

# julia> intlist
#  (2019-01-01, 2019-01-03)
#  (2019-01-01, 2019-01-06)
#  (2019-01-01, 2019-01-10)
#  (2019-01-03, 2019-01-06)
#  (2019-01-03, 2019-01-10)
#  (2019-01-06, 2019-01-10)

A = InsarTimeseries.build_A_matrix(geolist, intlist)
@test A == [1   0  0 ;
            0   1  0 ;
            0   0  1 ;
           -1   1  0 ;
           -1   0  1 ;
            0  -1  1 ]

B = InsarTimeseries.build_B_matrix(geolist, intlist)
@test B ≈ [2.0  0.0  0.0;
           2.0  3.0  0.0;
           2.0  3.0  4.0;
           0.0  3.0  0.0;
           0.0  3.0  4.0;
           0.0  0.0  4.0]

# TODO: make the const. vel into a function to test it
