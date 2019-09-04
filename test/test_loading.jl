using Test

demfile = joinpath(@__DIR__, "elevation.dem")

dem = InsarTimeseries.load(demfile)
rows, cols = (5, 6)

@test dem == reshape(collect(Int16, 1:30), rows, cols)
#  1   6  11  16  21  26
#  2   7  12  17  22  27
#  3   8  13  18  23  28
#  4   9  14  19  24  29
#  5  10  15  20  25  30

@test InsarTimeseries.load(demfile, (3, 4)) == dem[3, 4] == 18

@test_throws DomainError InsarTimeseries.load(demfile, (30, 40))
