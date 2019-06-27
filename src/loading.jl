import Glob
using PyCall


const SENTINEL_EXTS = [".geo", ".cc", ".int", ".amp", ".unw", ".unwflat"]
const COMPLEX_EXTS = [".int", ".slc", ".geo", ".cc", ".unw", ".unwflat", ".mlc", ".grd"]
const REAL_EXTS = [".amp", ".cor", ".mlc", ".grd"]
const ELEVATION_EXTS = [".dem", ".hgt"]
# These file types are not simple complex matrices: see load_stacked_img for detail
# .unwflat are same as .unw, but with a linear ramp removed
const STACKED_FILES = [".cc", ".unw", ".unwflat"]
const IMAGE_EXTS = [".png", ".tif", ".tiff", ".jpg"]

const LOAD_IN_PYTHON = vcat(ELEVATION_EXTS, IMAGE_EXTS, [".rsc", ".geojson", ".npy"])

const sario = PyNULL()

function __init__()
	copy!(sario, pyimport("apertools.sario"))
end

"""Examines file type for real/complex and runs appropriate load

Raises:
	ValueError: if sentinel files loaded without a .rsc file in same path
		to give the file width
"""
function load(filename::String; rsc_file::Union{String, Nothing}=nothing)
    ext = get_file_ext(filename)

	# For now, just pass through unimplemented extensions to Python
	if ext in LOAD_IN_PYTHON
		return sario.load(filename)
	end


    # Sentinel files should have .rsc file: check for dem.rsc, or elevation.rsc
    rsc_data = nothing
	if !isnothing(rsc_file)
		rsc_data = sario.load(rsc_file)::Dict{String, Any}
	end

    if ext in SENTINEL_EXTS
		if isnothing(rsc_file)
			rsc_file = sario.find_rsc_file(filename)
		end
        rsc_data = sario.load(rsc_file)
	end

    if ext in STACKED_FILES
        return load_stacked_img(filename, rsc_data)
	end
    # having rsc_data implies that this is not a UAVSAR file, so is complex
	# TODO: haven"t transferred over UAVSAR functions, so no load_real yet
	return load_complex(filename, rsc_data=rsc_data)

end

# # Make a shorter alias for load_file
# load = load_file

"""Extracts the file extension, including the "." (e.g.: .slc)"""
get_file_ext(filename::String) = splitext(filename)[end]


function read_stack(file_ext::String, rows::Int, cols::Int)
	allfiles = Glob.glob(file_ext, ".")
	buffer = Array{Float32, 2}(undef, (rows, cols))
    stack = Array{Float32, 3}(undef, (rows, cols, length(allfiles)))
    for (idx, f) in enumerate(allfiles)
 		read!(f, buffer)
 	   	stack[:, :, idx] = buffer
    end
    return stack
end

