import Glob
using HDF5

const SENTINEL_EXTS = [".geo", ".cc", ".int", ".amp", ".unw", ".unwflat"]
const COMPLEX_EXTS = [".int", ".slc", ".geo", ".cc", ".unw", ".unwflat", ".mlc", ".grd"]
const REAL_EXTS = [".amp", ".cor", ".mlc", ".grd"]
const ELEVATION_EXTS = [".dem", ".hgt"]
# These file types are not simple complex matrices: see load_stacked_img for detail
# .unwflat are same as .unw, but with a linear ramp removed
const STACKED_FILES = [".cc", ".unw", ".unwflat"]
const IMAGE_EXTS = [".png", ".tif", ".tiff", ".jpg"]

const LOAD_IN_PYTHON = vcat(IMAGE_EXTS, [".rsc", ".geojson", ".npy"])


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

    if ext in ELEVATION_EXTS
        return load_elevation(filename)
    end

    # Sentinel files should have .rsc file: check for dem.rsc, or elevation.rsc
    rsc_data = _get_rsc_data(filename, rsc_file)

    if ext in STACKED_FILES
        return load_stacked_img(filename, rsc_data)
    end
    # having rsc_data implies that this is not a UAVSAR file, so is complex
    # TODO: haven"t transferred over UAVSAR functions, so no load_real yet
    load_complex(filename, rsc_data)

end

"""Load one element of a file on disk (avoid reading in all of huge file"""
# TODO: Load a chunk of a file now?
function load(filename::String, row_col::Tuple{Int, Int}; rsc_file::Union{String, Nothing}=nothing)
    data_type = _get_data_type(filename)

    rsc_data = _get_rsc_data(filename, rsc_file)
    num_rows, num_cols = rsc_data["file_length"], rsc_data["width"]

    row, col = row_col
    if row < 1 || col < 1 || row > num_rows || col > num_cols
        throw(DomainError((row, col), " out of bounds for $filename of size ($num_rows, $num_cols)"))
    end
    seek_pos = _get_seek_position(row, col, num_cols, data_type)
    
    open(filename) do f
        seek(f, seek_pos)
        # This read syntax loads single `data_type`
        return read(f, data_type)
    end
end

"""For single element reading in binary files, seek to the right row, col"""
_get_seek_position(row, col, num_cols, data_type) = sizeof(data_type) * ((col - 1) + (num_cols * (row - 1)) )


function _get_rsc_data(filename, rsc_file)
    ext = get_file_ext(filename)

    rsc_data = nothing
    if !isnothing(rsc_file)
        rsc_data = sario.load(rsc_file)
    elseif ext in SENTINEL_EXTS || ext in ELEVATION_EXTS
        rsc_file = sario.find_rsc_file(filename)
        rsc_data = sario.load(rsc_file)
    end

    if !isnothing(rsc_data)
        rsc_data = convert(Dict{String, Any}, rsc_data)
    end
    return rsc_data
end

function _get_data_type(filename)
    ext = get_file_ext(filename)
    if ext in ELEVATION_EXTS
        return Int16
    elseif ext in STACKED_FILES
        return Float32
    else
        return ComplexF32
    end
end


"""Loads a digital elevation map from either .hgt file or .dem

.hgt is the NASA SRTM files given. Documentation on format here:
https://dds.cr.usgs.gov/srtm/version2_1/Documentation/SRTM_Topo.pdf
Key point: Big-endian 2 byte (16-bit) integers

.dem is format used by Zebker geo-coded and ROI-PAC SAR software
Only difference is data is stored little-endian (like other SAR data)

Note on both formats: gaps in coverage are given by INT_MIN -32768,
so either manually set data(data == np.min(data)) = 0,
or something like 
    data = clamp(data, -10000, Inf)
"""
function load_elevation(filename)
    ext = get_file_ext(filename)
    data_type = Int16 

    if ext == ".dem"
        rsc_file = sario.find_rsc_file(filename)
        dem_rsc = sario.load(rsc_file)
        rows, cols = (dem_rsc["file_length"], dem_rsc["width"])
        data = Array{data_type, 2}(undef, (cols, rows))

        read!(filename, data)

        # # TODO: Verify that the min real value will be above -1000
        # min_valid = -10000
        # # Set NaN values to 0
        # @. data[data < min_valid] = 0
        return transpose(data)
    else
        # swap_bytes = (ext == ".hgt")
        throw("$ext not implemented")
    end
end

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


function save_hdf5_stack(h5file::String, dset_name::String, stack; overwrite::Bool=true, do_permute=true)
    if overwrite
        mode = "w"
    else
        mode = "cw"  # Append mode
    end
    h5open(h5file, mode) do f 
        if do_permute
            write(f, dset_name, permutedims(stack, (2, 1, 3)))
        else
            write(f, dset_name, stack)
        end
    end
end

function save_deformation(h5file, 
                          deformation, 
                          geolist::Array{Date, 1}, 
                          dem_rsc; 
                          dset_name=STACK_DSET, 
                          do_permute=true,
                          unw_stack_file="unw_stack.h5")
    save_hdf5_stack(h5file, dset_name, deformation, do_permute=do_permute)
    h5open(h5file, "cw") do f
        geo_strings = Dates.format.(geolist, DATE_FMT)
        write(f, GEOLIST_DSET, geo_strings)
        # Note: using the sario python version since for some reason (as of 7.21.2019)
        # saving a JSON.json(dem) cant be loaded at all in h5py. weird string stuff
        # TODO: save reference that was used at the time
    end
    sario.save_dem_to_h5(h5file, dem_rsc, dset_name=DEM_RSC_DSET, overwrite=true)
end

function save_reference(h5file, unw_stack_file, dset_name, stack_flat_shifted_dset)
    # Read ref. unfo from unw file, save what was used to deformation result file
    reference = h5readattr(unw_stack_file, stack_flat_shifted_dset)[REFERENCE_ATTR]
    h5writeattr(h5file, dset_name, Dict(REFERENCE_ATTR => reference))
    reference_station = get(h5readattr(unw_stack_file, stack_flat_shifted_dset), REFERENCE_STATION_ATTR, "")
    h5writeattr(h5file, dset_name, Dict(REFERENCE_STATION_ATTR => reference_station))
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
