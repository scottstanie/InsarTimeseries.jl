using Dates

# Test data is 3 days after a starting reference on Jan 1:
# The actual phi phase is [0, 1, 2, 3]
 
geolist = [Date(2019,1,1), Date(2019,1,3), Date(2019,1,6), Date(2019,1,10)]
timediffs = InsarTimeseries.day_diffs(geolist)
@test timediffs == [2, 3, 4]

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

truth_phi = [0.0, 1.0, 2.0, 3.0]
truth_v = [1/2, 1/3, 1/4]
dphis = [1.0, 2.0, 3.0, 
         1.0, 2.0,
         1.0]
# Since we add no noise here, we should perfectly recover the phi phase
@test isapprox(A \ dphis, truth_phi[2:end], rtol=1e-7)

v_hat = B \ dphis
@test isapprox(v_hat, truth_v, rtol=1e-7)

phi_hat = InsarTimeseries.integrate_velocities(v_hat, timediffs)
@test isapprox(phi_hat, truth_phi, rtol=1e-7)

dpr = dphis .+ rand(eltype(dphis), size(dphis))
@test isapprox(A \ dpr,
               InsarTimeseries.integrate_velocities(B \ dpr, timediffs)[2:end],
               rtol=1e-7)
# TODO: make the const. vel into a function to test it
