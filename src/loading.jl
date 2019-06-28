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
		rsc_data = sario.load(rsc_file)
	end

    if ext in SENTINEL_EXTS
		if isnothing(rsc_file)
			rsc_file = sario.find_rsc_file(filename)
		end
		rsc_data = sario.load(rsc_file)
	end

	if !isnothing(rsc_data)
		rsc_data = convert(Dict{String, Any}, rsc_data)
	end

    if ext in STACKED_FILES
        return load_stacked_img(filename, rsc_data)
	end
    # having rsc_data implies that this is not a UAVSAR file, so is complex
	# TODO: haven"t transferred over UAVSAR functions, so no load_real yet
	load_complex(filename, rsc_data)

end

# # Make a shorter alias for load_file
# load = load_file

"""Extracts the file extension, including the "." (e.g.: .slc)"""
get_file_ext(filename::String) = splitext(filename)[end]

function load_complex(filename::String, rsc_data::Dict{String, Any})
	rows = rsc_data["file_length"]
	cols = rsc_data["width"]
	# Note: must be saved in c/row-major order, so loading needs a transpose
	out = Array{ComplexF32, 2}(undef, (cols, rows))

	read!(filename, out)
	return permutedims(out)
end


"""load_stacked_img is for weirdly formatted images:

Format is two stacked matrices:
	[[first], [second]] where the first "cols" number of floats
	are the first matrix, next "cols" are second, etc.
For .unw height files, the first is amplitude, second is phase (unwrapped)
For .cc correlation files, first is amp, second is correlation (0 to 1)
"""
function load_stacked_img(filename::String, rsc_data::Dict{String, Any})
	rows = rsc_data["file_length"]
	cols = rsc_data["width"]
	# Note: must be saved in c/row-major order, so loading needs a transpose
	# out_left usually is amplitude data
	# Usually we are interested in out_right
	#
	# First make container for all of data
	out = Array{Float32, 2}(undef, (2cols, rows))
	read!(filename, out)

	# TODO: port over rest of code for handling amp (if we care about that)
	# out_left = out[1:cols, :]
	return permutedims(out[cols+1:end, :])
end


# TODO: probably a better way to do this.. but can't figure out how to
# without allocating so many arrays that it's slow as Python
function load_stack(; file_list::Union{Array{String}, Nothing}=nothing, 
					directory::Union{String, Nothing}=nothing,
					file_ext::Union{String, Nothing}=nothing)
	if isnothing(file_list)
		file_list = sort(Glob.glob("*$file_ext", directory))
	end


	# rsc_data = sario.load(sario.find_rsc_file(basepath=directory))
	test_arr = load(file_list[1])
	rows, cols = size(test_arr)
	T = eltype(test_arr)

	stack, buffer = _return_array(T, rows, cols, length(file_list))
    for (idx, f) in enumerate(file_list)
 		read!(f, buffer)
 	   	stack[:, :, idx] = buffer
    end

	return _permute(stack, cols)
	# throw(ArgumentError(file_ext, "cant make stack from $file_ext"))
	# end
end

# Using multiple dispatch to avoid if statements for 2 types of images
function _return_array(T::Type{ComplexF32}, rows::Int, cols::Int, file_len::Int)
	stack = Array{T, 3}(undef, (cols, rows, file_len))
	buffer = Array{T, 2}(undef, (cols, rows))
	stack, buffer
end

function _return_array(T::Type{Float32}, rows::Int, cols::Int, file_len::Int)
	stack = Array{T, 3}(undef, (2cols, rows, file_len))
	buffer = Array{T, 2}(undef, (2cols, rows))
	stack, buffer
end

# For normal .int, .geo complex, just transpose stack
function _permute(stack::Array{ComplexF32, 3}, cols::Int)
	return permutedims(stack, (2, 1, 3))
end

# For normal weird stacked types, pick just the right half
function _permute(stack::Array{Float32, 3}, cols::Int)
	return permutedims(stack[cols+1:end, :, :], (2, 1, 3))
end



# function _load_stack_complex(file_list::Array{String}, rows::Int, cols::Int)
# 	stack = Array{ComplexF32, 3}(undef, (cols, rows, length(file_list)))
# 	buffer = Array{ComplexF32, 2}(undef, (cols, rows))
# 
#     for (idx, f) in enumerate(file_list)
#  		read!(f, buffer)
#  	   	stack[:, :, idx] = buffer
#     end
# 
# 	return permutedims(stack, (2, 1, 3))
# end
# 
# 
# function _load_stack_stacked(file_list::Array{String}, rows::Int, cols::Int)
# 	stack = Array{Float32, 3}(undef, (2cols, rows, length(file_list)))
# 	buffer = Array{Float32, 2}(undef, (2cols, rows))
# 
#     for (idx, f) in enumerate(file_list)
#  		read!(f, buffer)
#  	   	stack[:, :, idx] = buffer
#     end
# 	return permutedims(stack[cols+1:end, :, :], (2, 1, 3))
# end


# function save(filename::String, arr::Array{T, } ; kwargs...)
# sario.save(filename, arr, kwargs...)
# end
