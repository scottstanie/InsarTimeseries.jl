using Distributed: pmap, workers, WorkerPool

function prepare_stacks(igram_path; overwrite=false, ref_row=nothing,
                        ref_col=nothing, ref_station=nothing, 
                        geo_path=nothing, window=5, zero_masked=true)
    igram_path = abspath(igram_path)

    unw_stack_file = abspath(joinpath(igram_path, UNW_FILENAME))
    cc_stack_file = abspath(joinpath(igram_path, CC_FILENAME))
    mask_stack_file = abspath(joinpath(igram_path, MASK_FILENAME))

    # Step 1: create the masks.h5 file
    create_mask_stacks(igram_path, mask_filename=mask_stack_file,
                       overwrite=overwrite, geo_path=geo_path) 


    # Step 2: use these masks to zero bad ares to zero in .int and .cc files
    # Note: this may have been run already after making .int files, before unwrapping
    if zero_masked
        zero_masked_areas(igram_path, input_exts=[".int", ".cc", ".unw"])
    end

    create_hdf5_stack(cc_stack_file, ".cc", overwrite=overwrite, do_repack=true)

    deramp_unws(igram_path, input_ext=".unw", output_ext=".unwflat", overwrite=overwrite)

    create_hdf5_stack(unw_stack_file, ".unwflat", dset_name=STACK_FLAT_DSET, overwrite=overwrite)

    # TODO: do we need this after I have "NaN"s in the deramp function?
    # zero_masked_areas(igram_path, input_exts=[".unwflat"])

    if isnothing(ref_station) && isnothing(ref_row) && isnothing(ref_col)
    # TODO: auto pick, redo that logic here
    # ref_row, ref_col, ref_station = find_reference_location(
    #     unw_stack_file=unw_stack_file,
    #     cc_stack_file=cc_stack_file,
    #     mask_stack_file=mask_stack_file,
    # )
        ref_row, ref_col, ref_station = (50, 50, "")
    end

    shift_unw_file(unw_stack_file,
                   ref_row=ref_row,
                   ref_col=ref_col,
                   ref_station=ref_station,
                   window=window,
                   overwrite=overwrite)
end

# Note: 6 workers for the .geo.mask file creation seems fastest
# More and they thrash while reading
function create_mask_stacks(igram_path; mask_filename=nothing, geo_path=nothing, overwrite=false)
    if isnothing(mask_filename)
        mask_filename = abspath(joinpath(igram_path, MASK_FILENAME))
    end
    if isnothing(geo_path)
        geo_path = abspath(utils.get_parent_dir(igram_path))
    end

    row_looks, col_looks = utils.find_looks_taken(igram_path, geo_path=geo_path)
    dem_rsc = load(sario.find_rsc_file(directory=igram_path))

    loop_over_files(_get_geo_mask, geo_path, ".geo", ".geo.mask",
                    looks=(row_looks, col_looks), out_dir=igram_path,
                    max_procs=4, overwrite=overwrite)
                  
    # Now with bigger geo files read in parallel, 
    # write all to hdf5 file as stack
    save_masks(igram_path, geo_path, overwrite=overwrite)

    # Finall, add the aux. information
    dem_rsc = sario.load(sario.find_rsc_file(directory=igram_path))
    sario.save_dem_to_h5(mask_filename, dem_rsc, dset_name=DEM_RSC_DSET,
                         overwrite=overwrite)
    sario.save_geolist_to_h5(igram_path, mask_filename, overwrite=overwrite)
    sario.save_intlist_to_h5(igram_path, mask_filename, overwrite=overwrite)
end

_get_geo_mask(arr) = (arr .== 0)

"""Creates igram masks by taking the logical-or of the two .geo files

Assumes save_geo_masks already run"""
function save_masks(igram_path, geo_path; overwrite=false,
                    mask_filename=MASK_FILENAME,
                    geo_dset_name=GEO_MASK_DSET,
                    igram_dset_name=IGRAM_MASK_DSET)
    # TODO?: add checks for overwrite
    !sario.check_dset(mask_filename, geo_dset_name, overwrite) && return
    !sario.check_dset(mask_filename, igram_dset_name, overwrite) && return
    !sario.check_dset(mask_filename, GEO_MASK_SUM_DSET, overwrite) && return

    geo_mask_stack = load_stack(directory=igram_path,
                                file_ext=".geo.mask")


    int_date_list = sario.find_igrams(directory=igram_path)
    int_file_list = sario.find_igrams(directory=igram_path, parse=false)
    geo_date_list = sario.find_geos(directory=geo_path)

    nrows, ncols, _ = size(geo_mask_stack)
    int_mask_stack = Array{eltype(geo_mask_stack), 3}(undef, (nrows, ncols, length(int_file_list)))
    for (idx, (early, late)) in enumerate(int_date_list)
        early_idx = findfirst(isequal(early), geo_date_list)
        late_idx = findfirst(isequal(late), geo_date_list)
        early_mask = geo_mask_stack[:, :, early_idx]
        late_mask = geo_mask_stack[:, :, late_idx]

        out_filename = int_file_list[idx] * ".mask"

        new_mask = early_mask .| late_mask
        int_mask_stack[:, :, idx] .= new_mask
        # save(out_filename, new_mask)
    end
    save_hdf5_stack(mask_filename, igram_dset_name, convert(Array{Bool}, int_mask_stack), overwrite=false)
    save_hdf5_stack(mask_filename, geo_dset_name, convert(Array{Bool}, geo_mask_stack), overwrite=false)
            # write(f, dset_name, permutedims(stack, (2, 1, 3)))
    # Also create one image of the total masks
    # TODO: any way to have sum reduce the dims?
    h5write(mask_filename, GEO_MASK_SUM_DSET, 
            permutedims(sum(geo_mask_stack, dims=3)[:, :, 1]))
end

function _remaining_files(in_files::Array{String, 1}, out_files::Array{String, 1}, overwrite::Bool)
    if overwrite
        return in_files, out_files, collect(1:length(in_files))
    end
    in_todo = Array{String, 1}()
    out_todo = similar(in_todo)
    for (i, o) in zip(in_files, out_files)
        if !isfile(o)
            push!(in_todo, i)
            push!(out_todo, o)
        end
    end
    return in_todo, out_todo, indexin(in_todo, in_files)
end


# Step 2 functions: zeroing bad areas
function zero_masked_areas(directory; input_exts=[".int", ".cc"])
    for input_ext in input_exts
        in_files = find_files(input_ext, directory)
        out_files = in_files  # We are saving over the same file, just setting pixels to 0
        idxs = collect(1:length(in_files))

        # TODO: should we check for 0s in masked areas? maybe in the real function
        # in_todo, out_todo, idxs_todo = _remaining_files(in_files, out_files)

        wp = _get_workerpool(8)
        println("Starting zeroing areas for $input_ext")
        @time pmap((name_in, name_out, mask_idx) -> _load_and_run(zero_file, name_in, name_out, mask_idx,
                                                                  return_amp=true, do_permute=false),
                   wp, in_files, out_files, idxs)
        println("completed zeroing for $input_ext")
    end
end

# TODO: if i wanna make these true, default, i need to redo the _load_and_run arg passing
function zero_file(arr3d::AbstractArray{T, 3}, mask_idx; do_permute=false) where {T <: Number}
    mask = _read_mask(mask_idx, do_permute)
    data = @view arr3d[:, :, 2]
    data[mask] .= 0
    return arr3d
end

function zero_file(arr::AbstractArray{T, 2}, mask_idx; do_permute=false) where {T <: Number}
    mask = _read_mask(mask_idx, do_permute)
    arr[mask] .= 0
    return arr
end


function _read_mask(mask_idx, do_permute=false)
    m = Bool.(h5read(MASK_FILENAME, IGRAM_MASK_DSET, (:, :, mask_idx))[:, :, 1])
    return do_permute ? permutedims(m) : m
end


"""Remove a linear ramp from .unw files, save as .unwflat"""
function deramp_unws(directory="."; input_ext=".unw", output_ext=".unwflat", overwrite=true, do_permute=false)
    println("Writing deramped unws to $output_ext files")
    in_files = find_files(input_ext, directory)
    out_files = [replace(fname, input_ext => output_ext) for fname in in_files]

    # If we already have these files made, skip the function
    in_todo, out_todo, idxs_todo = _remaining_files(in_files, out_files, overwrite)

    println("$(length(out_todo)) files remain, overwrite = $overwrite")
    out_todo = overwrite ? in_todo : out_todo  # If we have overwrite=true, 
    if length(out_todo) == 0 && !overwrite
        println("skipping")
        return nothing
    end


    wp = _get_workerpool(8)
    println("Starting stack deramp ")
    pmap((name_in, name_out, mask_idx) -> _load_and_run(remove_ramp, name_in, name_out, mask_idx, do_permute=do_permute),
        wp, in_todo, out_todo, idxs_todo)
    println("Deramping stack complete")
end

sliceview(stack) = [view(stack, :, :, i) for i in 1:size(stack, 3)]


function remove_ramp(z, mask::AbstractArray{<:Number})
    # z_masked = copy(z)  # Do I really want a copy??
    z[mask] .= NaN
    return z - estimate_ramp(z)
end

function remove_ramp(z, mask_idx::Int; do_permute=false)
    mask = _read_mask(mask_idx, do_permute)
    return remove_ramp(z, mask)
end


"""Takes a 2D array an fits a linear plane to the data
    Ignores pixels that have nan or missing values

Since it is an order 1 surface, it will be 3 numbers, a, b, c from
     ax + by + c = z
"""
# TODO: Fix the mem-allocations here
function estimate_ramp(z::AbstractArray{<:AbstractFloat, 2})
    A = ones((length(z), 3))
    coeffs = Array{eltype(A), 1}(undef, (3,))
    z_fit = similar(z)

    zidxs = CartesianIndices(z)
    for idx in 1:size(A, 1)
        # if good_idxs[idx]
        # row, col is equiv to y, x, but subtract 1 to start at 0
        y, x = zidxs[idx].I .- 1
        A[idx, 1] = x
        A[idx, 2] = y
        # end
    end

    good_idxs = .~isnan.(z)
    coeffs .= A[vec(good_idxs), :] \ z[vec(good_idxs)]
    a, b, c = coeffs

    for idx in CartesianIndices(z_fit)
        # again subtract 1 to keep consistent with the fitting
        y, x = idx.I .- 1
        z_fit[idx] = a*x + b*y + c
    end

    return z_fit
end

"""Make stack as hdf5 file from a group of existing files"""
function create_hdf5_stack(filename::String,
                           file_ext::String;
                           dset_name=STACK_DSET,
                           directory=".",
                           overwrite=false,
                           do_repack=false)

    filename = isabspath(filename) ? filename : abspath(joinpath(directory, filename))
    println("Creating stack file $filename")

    get_file_ext(filename) in (".h5", ".hdf5") || ArgumentError("filename must end in .h5 or .hdf5")
    !sario.check_dset(filename, dset_name, overwrite) && return

    file_list = find_files(file_ext, directory)

    testf = sario.load(file_list[1])
    rows, cols = size(testf)
    # Note: since we're just loading/saving, don't permute on way in or way out
    shape = (cols, rows, length(file_list))
    _create_dset(filename, dset_name, shape, eltype(testf))
    arr_buf = zeros(eltype(testf), (cols, rows))
    mean_buf = zeros(eltype(testf), (cols, rows))

    h5open(filename, "r+") do hf
        dset = hf[dset_name]
        for (idx, f) in enumerate(file_list)
            arr_buf .= load(f, do_permute=false)
            dset[:, :, idx] = arr_buf
            mean_buf += arr_buf
        end
    end
    h5writeattr(filename, dset_name, Dict("filenames" => file_list))

    # Now save dem rsc as well
    dem_rsc = sario.load(sario.find_rsc_file(directory=directory))
    sario.save_dem_to_h5(filename, dem_rsc, dset_name=DEM_RSC_DSET, overwrite=overwrite)
    sario.save_geolist_to_h5(directory, filename, overwrite=overwrite)
    sario.save_intlist_to_h5(directory, filename, overwrite=overwrite)

    !sario.check_dset(filename, STACK_MEAN_DSET, overwrite) && return
    mean_buf ./= length(file_list)
    h5write(filename, STACK_MEAN_DSET, mean_buf)

    if do_repack
        repack(filename, dset_name)
    end

    return filename
end

function repack(filename, dset_name)
    # Now repack so that chunks are depth-wise (pixels for quick loading)
    nrows, ncols, nlayers = size(filename, dset_name)
    tmp_name = _make_tmp_name(filename)
    chunk_size = "$(nlayers)x10x10"
    println("Repacking $filename into size $chunk_size depth-wise chunks")
    cmd = `h5repack -v -f NONE -l $dset_name:CHUNK=$chunk_size $filename $tmp_name`
    println("Running:")
    println(cmd)
    @time run(cmd)
    run(`mv $tmp_name $filename`)
end
function _make_tmp_name(filename)
    fparts = splitpath(filename)
    return joinpath(fparts[1:end-1]..., "tmp_" * fparts[end])
end



"""Subtracts reference pixel group from each layer

window is the size around the reference pixel to 
average at each layer"""
function shift_unw_file(unw_stack_file::String; stack_flat_dset=nothing,
                        ref_row=nothing, ref_col=nothing, window=5, ref_station=nothing, 
                        overwrite::Bool=false, do_repack::Bool=true)
    """Runs a reference point shift on flattened stack of unw files stored in .h5
    
    do_repack will run h5repack to chunk"""
    if isnothing(stack_flat_dset)
        stack_flat_dset = STACK_FLAT_DSET
    end
    stack_flat_shifted_dset = stack_flat_dset * "_shifted"

    !sario.check_dset(unw_stack_file, stack_flat_shifted_dset, overwrite) && return

    if (isnothing(ref_row) || isnothing(ref_col))
        isnothing(ref_station) && throw(ArgumentError("Need ref_station if no ref_row/ref_col"))
        println("Using $ref_station as reference")
        rsc_data = Sario.load_dem_from_h5(unw_stack_file)
        ref_row, ref_col = gps.station_rowcol(station_name=ref_station, rsc_data=rsc_data)
    end
    println("Starting shift_stack: using ($ref_row, $ref_col) as reference")

    stack_size = size(unw_stack_file, stack_flat_dset)
    h5open(unw_stack_file, "r+") do f
        if !(stack_flat_dset in names(f))
            throw("Need $stack_flat_dset to be created in $unw_stack_file before shift stack can be run")
        end

        stack_in = f[stack_flat_dset]
        d_create(f,
            stack_flat_shifted_dset,
            datatype(Float32),
            dataspace(stack_size),
            # "chunk", (10, 10, size(stack_in, 3)),  # Seems better to repack?
            # Chunking here and writing by images is like 100x slower
        )
        stack_out = f[stack_flat_shifted_dset]

        # shift_stack(stack_in, stack_out, ref_row, ref_col, window=window)
        # Note: switching these so we don't have to permute dims upon loading HDF5Dset
        shift_stack(stack_in, stack_out, ref_col, ref_row, window=window)

    end

    if do_repack
        @time repack(unw_stack_file, stack_flat_shifted_dset)
    end

    h5writeattr(unw_stack_file, stack_flat_shifted_dset, Dict(REFERENCE_ATTR => [ref_row, ref_col]))
    # Make sure we don't write Nothing
    station_str = isnothing(ref_station) ? "" : ref_station
    h5writeattr(unw_stack_file, stack_flat_shifted_dset, Dict(REFERENCE_STATION_ATTR => station_str))

    println("Shifting stack complete")
end

"""
Note: window is size of the group around ref pixel to avg for reference.
    if window=1 or nothing, only the single pixel used to shift the group.
"""
# TODO: check later if worth doing one for HDF5Dataset, one for Array
function shift_stack(stack_in, stack_out, ref_row::Int, ref_col::Int; 
                     window::Int=5)
    half_win = div(window, 2)
    winsize = (2half_win + 1, 2half_win + 1)  # Make sure it is odd sized
    patch = Array{Float32, 2}(undef, winsize)

    nrows, ncols, nlayers = size(stack_in)
    layer = Array{Float32, 2}(undef, (nrows, ncols))
    layer_out = similar(layer)

    @show stack_in
    @show stack_out

    @inbounds for k = 1:nlayers
        layer .= view(stack_in, :, :, k)
        layer_out .=  _shift_layer(layer, patch, ref_row, ref_col, half_win)
        stack_out[:, :, k] = layer_out

        if k % 100 == 0
            println("Finished with $k layers")
        end
    end
    return stack_out
end
function shift_stack(stack_in, ref_row::Int, ref_col::Int; window::Int=5)
    stack_out = Array{eltype(stack_in), ndims(stack_in)}(undef, size(stack_in))
    return shift_stack(stack_in, stack_out, ref_row, ref_col; window=window)
end

function _shift_layer(layer, patch, ref_row, ref_col, half_win)
    patch .= layer[ref_row - half_win:ref_row + half_win, 
                   ref_col - half_win:ref_col + half_win]

    # Adding the `view` to eliminate extra singleton dimension from HDF5
    return view(layer .- mean(patch), :, :)
end


function create_mean_hdf5(h5file::String; dset_name::String=STACK_DSET)
    h5open(h5file, "cw") do f
        dset = f[dset_name]
        nrows, ncols, nlayers = size(dset)
        mean_buf = zeros(eltype(dset), nrows, ncols)

        println("Creating mean from stack:")
        @show size(dset)
        for idx in 1:nlayers
            # Adding the `[:,:]` to eliminate extra singleton dimension from HDF5
            mean_buf .+= dset[:, :, idx][:, :]
        end
        mean_buf ./= nlayers
        write(f, STACK_MEAN_DSET, mean_buf)
    end
end


function _create_dset(h5file, dset_name, shape, dtype)
    # TODO: add a chunking option for ones we access depth-wise vs layer-wise
    # dset = d_create(g, "B", datatype(Float64), dataspace(1000,100,10), "chunk", (100,100,1))
    h5open(h5file, "r+") do f
        d_create(f,
            dset_name,
            datatype(dtype),
            dataspace(shape),
        )
    end
end



"""Parallel way to run over many image files
Pass a function `f` to operate on each image 

    With 8 workers, on 100 files of 1000x1000

    @time InsarTimeseries.loop_over_files(add1, ".", ".int", "_flat.int")

    runs in 1 second, vs about 8 seconds serially
"""
function loop_over_files(f, directory::String, input_ext::String, output_ext::String;
                         looks=(1, 1), out_dir=nothing,
                         max_procs::Union{Nothing, Int}=nothing,
                         overwrite=false, do_permute=true)
    in_files = find_files(input_ext, directory)
    println("Looping over files: $in_files")
    out_files = [replace(fname, input_ext => output_ext) for fname in in_files]
    if !isnothing(out_dir)
        # Redo the path out if specified
        out_files = [joinpath(out_dir, splitpath(f)[end]) for f in out_files]
    end

    in_todo, out_todo, idxs_todo = _remaining_files(in_files, out_files, overwrite)
    println("$(length(out_todo)) files remain, overwrite = $overwrite")
    if length(out_todo) == 0 && !overwrite
        println("skipping")
        return nothing
    end

    wp = _get_workerpool(max_procs)
    pmap((name_in, name_out) -> _load_and_run(f, name_in, name_out, looks=looks, do_permute=do_permute),
         wp, in_todo, out_todo)
    # @sync @distributed for fname in in_files
    # end
end

"""To manually pass which files you want to run on and save to""" 
function loop_over_files(f, in_files::Array{String, 1}, out_files::Array{String, 1};
                         looks=(1, 1), max_procs::Union{Nothing, Int}=nothing, do_permute=true)
    wp = _get_workerpool(max_procs)
    pmap((name_in, name_out) -> _load_and_run(f, name_in, name_out, looks=looks, do_permute=do_permute),
         wp, in_files, out_files)
end

function _load_and_run(f, name_in, name_out, args...; looks=(1, 1), kwargs...)
    println("Input: $name_in, out: $name_out")
    input_arr = load(name_in; looks=looks, kwargs...)
    out_arr = f(input_arr, args...)
    save(name_out, out_arr; kwargs...)
end

function _get_workerpool(max_procs)
    if isnothing(max_procs)
        return WorkerPool(workers())
    else
        m = min(max_procs, length(workers()))
        return WorkerPool(workers()[1:m])
    end
end
