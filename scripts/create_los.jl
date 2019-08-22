import Glob
println("Loading:")
@time include("../src/InsarTimeseries.jl")

if ispath("../extra_files/")
    db_file = Glob.glob("../extra_files/*.db*")[1]
else
    db_file = Glob.glob("../*.db*")[1]
end

dem_rsc = InsarTimeseries.sario.load("dem.rsc")

@time InsarTimeseries.los_map(dem_rsc, db_file, "los_map.h5");
