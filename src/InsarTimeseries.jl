"""Module InsarTimeseries

Testing docstring
"""
module InsarTimeseries

using Dates

include("./loading.jl")
include("./timeseries.jl")


export load,
	# load_file, 
	read_stack, 
	get_file_ext,
	matrixA

end # module

