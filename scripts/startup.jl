@time using HDF5
@time import InsarTimeseries, MapImages, Sario
repostr = expanduser("~/repos/")
@time include(repostr*"/InsarTimeseries.jl/scripts/plotting.jl")
@time include(repostr*"/InsarTimeseries.jl/scripts/find_abs_shift.jl")
@time include(repostr*"/InsarTimeseries.jl/scripts/point_analysis.jl")
