"""Module InsarTimeseries

Testing docstring
"""
module InsarTimeseries

using Dates

# Include sario here so it is available to all
using PyCall
const sario = PyNULL()

function __init__()
    copy!(sario, pyimport("apertools.sario"))
end

include("./loading.jl")
include("./timeseries.jl")
include("./mask.jl")


export load,
    run_inversion,
    # load_file, 
    load_stack, 
    get_file_ext

end # module

