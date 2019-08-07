"""Module InsarTimeseries

Testing docstring
"""
module InsarTimeseries

# Include apertools Python modeules here to make available to all
using PyCall
const sario = PyNULL()
const gps = PyNULL()

function __init__()
    copy!(sario, pyimport("apertools.sario"))
    copy!(gps, pyimport("apertools.gps"))
end


# Type alias for commonly used compositite type
const Igram = Tuple{Date, Date}


include("./common.jl")
include("./loading.jl")
include("./prepare.jl")
include("./stackavg.jl")
include("./sbas.jl")
include("./timeseries.jl")
include("./mask.jl")
include("./projections.jl")
include("./los.jl")


# TODO: figure out which we care to export
# export load,
#     run_inversion,
#     ...

end # module

