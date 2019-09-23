"""Module InsarTimeseries
"""

__precompile__(true)

module InsarTimeseries

# Include apertools Python modules here to make available to all
using PyCall
const sario = PyNULL()
const gps = PyNULL()
const utils = PyNULL()


#TODO: functions needed to remove these python calls:
# prepare.jl: a bunch of `sario`
# gps.load_gps_los_data
# scripts/find_abs_shift.jl:    lon, lat = gps.station_lonlat(station_name)
# scripts/testl1.jl:  lon, lat = gps.station_lonlat(station_name)
# src/analysis.jl: dts, gps_los_data = InsarTimeseries.gps.load_gps_los_data(geo_path, station_name,
# src/prepare.jl: ref_row, ref_col = gps.station_rowcol(station_name=ref_station, rsc_data=rsc_data)
# src/prepare.jl: geo_path = abspath(utils.get_parent_dir(igram_path))
# src/prepare.jl: row_looks, col_looks = utils.find_looks_taken(igram_path, geo_path=geo_path)
# srx/prepare.jl: utils.get_parent_dir(igram_path))
# src/prepare.jl:    row_looks, col_looks = utils.find_looks_taken
function __init__()
    copy!(sario, pyimport("apertools.sario"))
    copy!(gps, pyimport("apertools.gps"))
    copy!(utils, pyimport("apertools.utils"))
end

using Sario

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

