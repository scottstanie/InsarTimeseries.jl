"""Module InsarTimeseries
"""

__precompile__(true)

module InsarTimeseries

# Include apertools Python modules here to make available to all
using PyCall
const sario = PyNULL()
const gps = PyNULL()
const utils = PyNULL()

function __init__()
    copy!(sario, pyimport("apertools.sario"))
    copy!(gps, pyimport("apertools.gps"))
    copy!(utils, pyimport("apertools.utils"))
end

using Sario

include("./common.jl")
include("./prepare.jl")
include("./stackavg.jl")
include("./optimize.jl")
include("./sbas.jl")
include("./timeseries.jl")
include("./mask.jl")
include("./projections.jl")
include("./los.jl")
include("./analysis.jl")


# TODO: figure out which we care to export
export STACK_DSET,
       STACK_DSET,
       STACK_MEAN_DSET,
       STACK_FLAT_DSET,
       STACK_FLAT_SHIFTED_DSET,
       REFERENCE_ATTR,
       REFERENCE_STATION_ATTR
        #     run_inversion,
    #     ...

end # module

