"""Module InsarTimeseries
"""

__precompile__(true)

module InsarTimeseries

using Sario
import MapImages

include("./common.jl")
include("./config.jl")
include("./prepare.jl")
include("./stackavg.jl")
include("./optimize.jl")
include("./sbas.jl")
include("./timeseries.jl")
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
