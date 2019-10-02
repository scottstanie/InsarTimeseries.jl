import Glob
println("Loading:")
try
    import InsarLOS
catch
    println("insar los import error")
end
try
    import Sario
catch
    println("Sario import error")
end

if ispath("../extra_files/")
    dbpath = Glob.glob("../extra_files/*.db*")[1]
else
    dbpath = Glob.glob("../*.db*")[1]
end


@time InsarLOS.create_los_map(directory=".", dbpath=dbpath, outfile="los_map.h5");
