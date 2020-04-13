import Dates
import Sario

dfmt = dateformat"yyyymmdd"

_to_str(dpair::Tuple{Date,Date}, ext) =
    Dates.format(dpair[1], dfmt) * "_" * Dates.format(dpair[2], dfmt) * ext


function avg_files(
    geolist::AbstractArray,
    intlist::AbstractArray,
    ext::String = ".unw",
    outdir::String = "avg_igrams/",
)
    mkpath(outdir)  # Like mkdir -p, wont error

    nrows, ncols = size(Sario.load("dem.rsc"))
    for d in geolist
        buf = zeros(Float32, (ncols, nrows))
        curints = intlist[d in intlist]
        for iname in curints
            cur = Sario.load(_to_str(iname, ext), do_permute = false)
            cur[isnan.(cur)] .= 0
            buf .+= cur
        end
        buf ./= length(curints)
        Sario.save(joinpath(outdir, Dates.format(d, dfmt) * ext), buf, do_permute = false)
        println("Saved $d")
    end
end

function avg_files(ext::String = ".unw", outdir::String = "avg_igrams/")
    geolist = Sario.load_geolist_from_h5("unw_stack.h5")  # TODO: just glob from dir
    intlist = Sario.load_intlist_from_h5("unw_stack.h5")
    return avg_files(geolist, intlist, ext, outdir)
end
