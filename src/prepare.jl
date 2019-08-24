using Distributed: pmap, myid, workers, WorkerPool

find_files(ext, directory=".") = sort(Glob.glob("*"*ext, directory))

"""Runs a reference point shift on flattened stack of unw files stored in .h5"""
function deramp_unws(directory="."; input_ext=".unw", output_ext=".unwflat", overwrite=true)
    println("Writing deramped unws to $output_ext files")
    in_files = find_files(input_ext, directory)
    out_files = [replace(fname, input_ext => output_ext) for fname in in_files]

    # If we already have these files made, skip the function
    if all(isfile(f) for f in out_files) && !overwrite
        println("Output Files, overwrite = $overwrite. Skipping")
        return nothing
    end

    println("Starting stack deramp ")
    loop_over_files((arr, mask) -> remove_ramp(arr, mask), in_files, out_files)
    println("Deramping stack complete")
end

function _create_dset(h5file, dset_name, shape, dtype)
    # TODO: add a chunking option for ones we access depth-wise vs layer-wise
    # dset = d_create(g, "B", datatype(Float64), dataspace(1000,100,10), "chunk", (100,100,1))
    h5open(h5file, "cw") do f
        d_create(f,
            dset_name,
            datatype(dtype),
            dataspace(shape),
        )
    end
end

# Note: 6 workers for the .geo.mask file creation seems fastest
# More and they thrash while reading
function create_mask_stacks(igram_path, mask_filename=nothing, geo_path=nothing, overwrite=false)
    if isnothing(mask_filename)
        mask_file = joinpath(igram_path, MASK_FILENAME)
    end
    if isnothing(geo_path)
        geo_path = utils.get_parent_dir(igram_path)
    end

    row_looks, col_looks = utils.find_looks_taken(igram_path, geo_path=geo_path)
    dem_rsc = load(sario.find_rsc_file(directory=igram_path))
    # TODO: make sure the masks geo/intlists are togeher...
    # Maybe just lump them with the unw file
    # sario.save_dem_to_h5(mask_file, dem_rsc, dset_name=DEM_RSC_DSET, overwrite=overwrite)
    # save_geolist_to_h5(igram_path, mask_file, overwrite=overwrite)
    # save_intlist_to_h5(igram_path, mask_file, overwrite=overwrite)

    # TODO: do the overwrite check
    loop_over_files(_get_geo_mask, geo_path, ".geo", ".geo.mask",
                    looks=(row_looks, col_looks))

    compute_int_masks(
        mask_file=mask_file,
        igram_path=igram_path,
        geo_path=geo_path,
        row_looks=row_looks,
        col_looks=col_looks,
        dem_rsc=dem_rsc,
        overwrite=overwrite,
    )
end

_get_geo_mask(arr) = (arr .== 0)

"""Creates igram masks by taking the logical-or of the two .geo files

Assumes save_geo_masks already run
"""
function compute_int_masks(
        igram_path=None,
        geo_path=None,
        dem_rsc=None,
        dset_name=IGRAM_MASK_DSET,
        overwrite=False,
)
    # TODO: add checks for overwrite

    int_date_list = sario.find_igrams(directory=igram_path)
    int_file_list = sario.find_igrams(directory=igram_path, parse=false)

    geo_date_list = sario.find_geos(directory=geo_path)


    geo_mask_stack = f[GEO_MASK_DSET]
    int_mask_dset = f[dset_name]
    # TODO: fix/pass to parallel
    for (idx, (early, late)) in enumerate(int_date_list)
        early_idx = geo_date_list.index(early)
        late_idx = geo_date_list.index(late)
        early_mask = geo_mask_stack[early_idx]
        late_mask = geo_mask_stack[late_idx]

        new_mask = early_mask .| late_mask
        int# _mask_sum .+= new_mask

    end
    # # Also create one image of the total masks
    # f[IGRAM_MASK_SUM_DSET] = sum(int_mask_dset, dims=0)
end


"""Subtracts reference pixel group from each layer

window is the size around the reference pixel to average at each layer
"""
function shift_unw_file(unw_stack_file::String; stack_flat_dset=nothing,
                        ref_row=nothing, ref_col=nothing, window=5, ref_station=nothing, 
                        overwrite=false)
    """Runs a reference point shift on flattened stack of unw files stored in .h5"""
    if isnothing(stack_flat_dset)
        stack_flat_dset = STACK_FLAT_DSET
    end
    stack_flat_shifted_dset = stack_flat_dset * "_shifted"

    if !sario.check_dset(unw_stack_file, stack_flat_shifted_dset, overwrite)
        return nothing
    end

    if (isnothing(ref_row) || isnothing(ref_col))
        println("Using $ref_station as reference")
        rsc_data = sario.load_dem_from_h5(unw_stack_file)
        ref_row, ref_col = gps.station_rowcol(station_name=ref_station, rsc_data=rsc_data)
    end
    println("Starting shift_stack: using ($ref_row, $ref_col) as reference")

    h5open(unw_stack_file, "cw") do f
        if !(stack_flat_dset in names(f))
            throw("Need $stack_flat_dset to be created in $unw_stack_file before shift stack can be run")
        end

        stack_in = f[stack_flat_dset]
        d_create(f,
            stack_flat_shifted_dset,
            datatype(Float32),
            dataspace(size(stack_in)),
        )
        stack_out = f[stack_flat_shifted_dset]

        # shift_stack(stack_in, stack_out, ref_row, ref_col, window=window)
        # Note: switching these so we don't have to permute dims upon loading HDF5Dset
        shift_stack(stack_in, stack_out, ref_col, ref_row, window=window)

    end

    h5writeattr(unw_stack_file, stack_flat_shifted_dset, Dict(REFERENCE_ATTR => [ref_row, ref_col]))
    # Make sure we don't write Nothing
    station_str = isnothing(ref_station) ? "" : ref_station
    h5writeattr(unw_stack_file, stack_flat_shifted_dset, Dict(REFERENCE_STATION_ATTR => station_str))

    println("Shifting stack complete")
end


"""
Note: window is size of the group around ref pixel to avg for reference.
    if window=1 or None, only the single pixel used to shift the group.
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

function _shift_layer(layer, patch, ref_row, ref_col, half_win)
    patch .= layer[ref_row - half_win:ref_row + half_win, 
                   ref_col - half_win:ref_col + half_win]

    # Adding the `view` to eliminate extra singleton dimension from HDF5
    return view(layer .- mean(patch), :, :)
end


function remove_ramp(z, mask; buf=nothing, A=nothing, coeffs=nothing, z_fit=nothing)
    if isnothing(buf)
        z_masked = copy(z)
    else
        z_masked = buf
        z_masked .= z
    end

    z_masked[mask] .= NaN
    return z - estimate_ramp(z_masked, A, coeffs, z_fit)
end


"""Takes a 2D array an fits a linear plane to the data
    Ignores pixels that have nan or missing values

Since it is an order 1 surface, it will be 3 numbers, a, b, c from
     ax + by + c = z
"""
# TODO: Fix the mem-allocations here
function estimate_ramp(z::Array{<:AbstractFloat, 2}, A=nothing, coeffs=nothing, z_fit=nothing)
    good_idxs = .~isnan.(z)

    if isnothing(A)
        A = ones((length(z), 3))
    end
    if isnothing(coeffs)
        coeffs = Array{eltype(A), 1}(undef, (3,))
    end
    if isnothing(z_fit)
        z_fit = similar(z)
    end
    Aidx = 1

    for idx in CartesianIndices(z)
        if good_idxs[idx]
            # row, col is equiv to y, x, but subtract 1 to start at 0
            y, x = idx.I .- 1
            A[Aidx, 1] = x
            A[Aidx, 2] = y
            Aidx += 1
        end
    end
    coeffs .= A[good_idxs, :] \ z[good_idxs]
    a, b, c = coeffs

    for idx in CartesianIndices(z_fit)
        # again subtract 1 to keep consistent with the fitting
        y, x = idx.I .- 1
        z_fit[idx] = a*x + b*y + c
    end

    return z_fit
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


"""Parallel way to run over many image files
Pass a function `f` to operate on each image 

    With 8 workers, on 100 files of 1000x1000

    @time InsarTimeseries.loop_over_files(add1, ".", ".int", "_flat.int")

    runs in 1 second, vs about 8 seconds serially
"""
function loop_over_files(f, directory::String, input_ext::String, output_ext::String;
                         looks=(1, 1), max_procs::Union{Nothing, Int}=nothing)
    in_files = find_files(input_ext, directory)
    println("Looping over files: $in_files")
    out_files = [replace(fname, input_ext => output_ext) for fname in in_files]

    wp = _get_workerpool(max_procs)
    pmap((name_in, name_out) -> _load_and_run(f, name_in, name_out, looks=looks),
         wp, in_files, out_files)
    # @sync @distributed for fname in in_files
    # end
end

"""To manually pass which files you want to run on and save to""" 
function loop_over_files(f, in_files::Array{String, 1}, out_files::Array{String, 1};
                         looks=(1, 1), max_procs::Union{Nothing, Int}=nothing)
    wp = _get_workerpool(max_procs)
    pmap((name_in, name_out) -> _load_and_run(f, name_in, name_out, looks=looks),
         wp, in_files, out_files)
end

function _load_and_run(f, name_in, name_out; looks=(1, 1))
    println("On proc id $(myid()): Input: $name_in, out: $name_out")
    input_arr = load(name_in, looks=looks)
    out_arr = f(input_arr)
    save(name_out, out_arr)
end
function _get_workerpool(max_procs)
    if isnothing(max_procs)
        return WorkerPool(workers())
    else
        return WorkerPool(workers()[1:max_procs])
    end
end
