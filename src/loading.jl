import Glob
using HDF5

const DATE_FMT = "yyyymmdd"
const SENTINEL_EXTS = [".geo", ".cc", ".int", ".amp", ".unw", ".unwflat"]
const COMPLEX_EXTS = [".int", ".slc", ".geo", ".cc", ".unw", ".unwflat", ".mlc", ".grd"]
const REAL_EXTS = [".amp", ".cor", ".mlc", ".grd"]
const ELEVATION_EXTS = [".dem", ".hgt"]
# These file types are not simple complex matrices: see load_stacked_img for detail
# .unwflat are same as .unw, but with a linear ramp removed
const STACKED_FILES = [".cc", ".unw", ".unwflat"]
const IMAGE_EXTS = [".png", ".tif", ".tiff", ".jpg"]

const LOAD_IN_PYTHON = vcat(ELEVATION_EXTS, IMAGE_EXTS, [".rsc", ".geojson", ".npy"])

# dataset names for general 3D stacks
const STACK_DSET = "stack"
const STACK_MEAN_DSET = "mean_stack"
const STACK_FLAT_DSET = "deramped_stack"

# Mask file datasets
const GEO_MASK_DSET = "geo"
const IGRAM_MASK_DSET = "igram"
const IGRAM_MASK_SUM_DSET = "igram_sum"

const DEM_RSC_DSET = "dem_rsc"
const GEOLIST_DSET = "geo_dates"
const INTLIST_DSET = "int_dates"


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

function load_geolist_from_h5(h5file::String)
    h5open(h5file) do f
        geo_strings = read(f, GEOLIST_DSET)
        return sario.parse_geolist_strings(geo_strings)
    end
end


function load_intlist_from_h5(h5file)
    h5open(h5file) do f
        int_strings = read(f, INTLIST_DSET)
        # Note transpose, since it's stored as a N x 2 array
        # (which is flipped to 2 x N for julia
        return sario.parse_intlist_strings(int_strings')
    end
end


function save_hdf5_stack(h5file::String, dset_name::String, stack; overwrite::Bool=true)
    if overwrite
        mode = "w"
    else
        mode = "cw"  # Append mode
    end
    h5open(h5file, mode) do f 
        write(f, dset_name, permutedims(stack, (2, 1, 3)))
    end
end

function save_deformation(h5file, deformation, geolist::Array{Date, 1}, dem_rsc; dset_name=STACK_DSET)
    save_hdf5_stack(h5file, dset_name, deformation)
    h5open(h5file, "cw") do f
        geo_strings = Dates.format.(geolist, DATE_FMT)
        write(f, GEOLIST_DSET, geo_strings)
        # Note: using the sario python version since for some reason (as of 7.21.2019)
        # saving a JSON.json(dem) cant be loaded at all in h5py. weird string stuff
        sario.save_dem_to_h5(h5file, dem_rsc, dset_name=DEM_RSC_DSET, overwrite=true)
    end
end

"""Wrapper around h5read to account for the transpose
necessary when reading from python-written stacks

"""
function load_hdf5_stack(h5file::String, dset_name::String)
    h5open(h5file) do f
        return permutedims(read(f[dset_name]), (2, 1, 3))
    end
end

# If loading only certain layers, don't read all into memory
function load_hdf5_stack(h5file::String, dset_name::String, valid_layer_idxs)
    h5open(h5file) do f
        dset = f[dset_name]
        nrows, ncols, _ = size(dset)
        out = Array{eltype(dset), ndims(dset)}(undef, (nrows, ncols, length(valid_layer_idxs)))
        for (tot_idx, v_idx) in enumerate(valid_layer_idxs)
            out[:, :, tot_idx] = dset[:, :, v_idx]
        end
        return permutedims(out, (2, 1, 3))
    end
end

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
#     stack = Array{ComplexF32, 3}(undef, (cols, rows, length(file_list)))
#     buffer = Array{ComplexF32, 2}(undef, (cols, rows))
# 
#     for (idx, f) in enumerate(file_list)
#          read!(f, buffer)
#             stack[:, :, idx] = buffer
#     end
# 
#     return permutedims(stack, (2, 1, 3))
# end
# 
# 
# function _load_stack_stacked(file_list::Array{String}, rows::Int, cols::Int)
#     stack = Array{Float32, 3}(undef, (2cols, rows, length(file_list)))
#     buffer = Array{Float32, 2}(undef, (2cols, rows))
# 
#     for (idx, f) in enumerate(file_list)
#          read!(f, buffer)
#             stack[:, :, idx] = buffer
#     end
#     return permutedims(stack[cols+1:end, :, :], (2, 1, 3))
# end


# function save(filename::String, arr::Array{T, } ; kwargs...)
# sario.save(filename, arr, kwargs...)
# end
