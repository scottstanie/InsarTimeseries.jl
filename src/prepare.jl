"""Subtracts reference pixel group from each layer

window is the size around the reference pixel to average at each layer
"""
function shift_unw_file(unw_stack_file::String; ref_row=nothing, ref_col=nothing, window=3, ref_station=nothing, overwrite=false)
    """Runs a reference point shift on flattened stack of unw files stored in .h5"""
    if !sario.check_dset(unw_stack_file, STACK_FLAT_SHIFTED_DSET, overwrite)
        return nothing
    end

    if (isnothing(ref_row) || isnothing(ref_col))
        println("Using $ref_station as reference")
        rsc_data = sario.load_dem_from_h5(unw_stack_file)
        ref_row, ref_col = gps.station_rowcol(station_name=ref_station, rsc_data=rsc_data)
    end
    println("Starting shift_stack: using ($ref_row, $ref_col) as reference")

    h5open(unw_stack_file, "cw") do f
        if !(STACK_FLAT_DSET in names(f))
            throw("Need $STACK_FLAT_DSET to be created in $unw_stack_file before shift stack can be run")
        end

        stack_in = f[STACK_FLAT_DSET]
        d_create(f,
            STACK_FLAT_SHIFTED_DSET,
            datatype(Float32),
            dataspace(size(stack_in)),
        )
        stack_out = f[STACK_FLAT_SHIFTED_DSET]

        # shift_stack(stack_in, stack_out, ref_row, ref_col, window=window)
        # Note: switching these so we don't have to permute dims upon loading HDF5Dset
        shift_stack(stack_in, stack_out, ref_col, ref_row, window=window)

    end

    h5writeattr(unw_stack_file, STACK_FLAT_SHIFTED_DSET, Dict(REFERENCE_ATTR => [ref_row, ref_col]))
    # Make sure we don't write Nothing
    station_str = isnothing(ref_station) ? "" : ref_station
    h5writeattr(unw_stack_file, STACK_FLAT_SHIFTED_DSET, Dict(REFERENCE_STATION_ATTR => station_str))

    println("Shifting stack complete")
end


"""
Note: window is size of the group around ref pixel to avg for reference.
    if window=1 or None, only the single pixel used to shift the group.
"""

# TODO: check later if worth doing one for HDF5Dataset, one for Array
function shift_stack(stack_in, stack_out, ref_row::Int, ref_col::Int; 
                     window::Int=3, chunk_layers=100)
    half_win = div(window, 2)
    winsize = (2half_win + 1, 2half_win + 1)  # Make sure it is odd sized
    patch = Array{Float32, 2}(undef, winsize)

    @show stack_in
    @show stack_out

    # Read chunks in at a time to limit HDF5 reading latency
    nrows, ncols, nlayers = size(stack_in)
    chunk =  Array{Float32, 3}(undef, (nrows, ncols, chunk_layers))

    k = 1
    while k < size(stack_in, 3)
        endk = min(k + chunk_layers - 1, lastindex(stack_in, 3))
        println("Shiffting layers $k : $endk")
        chunk[:, :, k:endk] .= stack_in[:, :, k:endk]

        for jdx in 1:(endk - k)
            layer = view(chunk, :, :, jdx)
        
            patch .= layer[ref_row - half_win:ref_row + half_win, 
                           ref_col - half_win:ref_col + half_win]

            # Adding the `view` to eliminate extra singleton dimension from HDF5
            stack_out[:, :, k + jdx - 1] = view(layer .- mean(patch), :, :)
        end
        k += chunk_layers
    end
    return stack_out
end

function deramp_unw_file(unw_stack_file::String; order=2, overwrite=false)
    """Runs a reference point shift on flattened stack of unw files stored in .h5"""
    if order == 1
        flat_dset = STACK_FLAT_DSET1
    elseif order == 2
        flat_dset = STACK_FLAT_DSET2
    end

    if !sario.check_dset(unw_stack_file, flat_dset, overwrite)
        return nothing
    end

    println("Starting stack deramp : using order $order")
    fmask = h5open("masks.h5", "r")
    h5open(unw_stack_file, "cw") do f
        if !(STACK_DSET in names(f))
            throw("Need $STACK_DSET to be created in $unw_stack_file before deramp stack can be run")
        end

        stack_in = f[STACK_DSET]
        d_create(f,
            flat_dset,
            datatype(Float32),
            dataspace(size(stack_in)),
        )
        stack_out = f[flat_dset]
        mask_dset = fmask[IGRAM_MASK_DSET]

        for k = 1:size(stack_in, 3)
            layer = view(stack_in, :, :, k)
            mask = Bool.(view(mask_dset, :, :, k))
            # @show size(layer), size(mask), size(stack_out[:,:,k])
            stack_out[:, :, k] = remove_ramp(layer, order, mask)
        end

    end
    close(fmask)

    h5writeattr(unw_stack_file, STACK_FLAT_SHIFTED_DSET, Dict("order" => order))

    println("Deramping stack complete")
end

function remove_ramp(z, order, mask)
    z_masked = copy(z)
    z_masked[mask] .= NaN
    if order == 1
        return z - estimate_ramp1(z_masked)
    elseif order == 2
        return z - estimate_ramp2(z_masked)
    else
        println("WARNING: Order $order not supported. Running order 1 ramp removal")
        return z - estimate_ramp1(z_masked)
    end

end


"""Takes a 2D array an fits a linear plane to the data
    Ignores pixels that have nan or missing values

For order = 1, it will be 3 numbers, a, b, c from
     ax + by + c = z
"""
function estimate_ramp1(z::Array{<:AbstractFloat, 2})
    good_idxs = .~isnan.(z)

    A = ones((sum(good_idxs), 3))
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
    coeffs = A \ z[good_idxs]
    a, b, c = coeffs

    z_fit = similar(z)
    for idx in CartesianIndices(z_fit)
        # again subtract 1 to keep consistent with the fitting
        y, x = idx.I .- 1
        z_fit[idx] = a*x + b*y + c
    end

    return z_fit
end

"""
For order = 2, it will be 6:
    a*x + b*y + c*x*y + d*x^2 + ey^2 + f
"""
function estimate_ramp2(z::Array{<:AbstractFloat, 2})
    good_idxs = .~isnan.(z)

    A = ones((sum(good_idxs), 6))
    Aidx = 1

    for idx in CartesianIndices(z)
        if good_idxs[idx]
            # row, col is equiv to y, x, but subtract 1 to start at 0
            y, x = idx.I .- 1
            A[Aidx, 1:5] = [x, y, x*y, x^2, y^2]
            Aidx += 1
        end
    end
    coeffs = A \ z[good_idxs]
    a, b, c, d, e, f = coeffs

    z_fit = similar(z)
    for idx in CartesianIndices(z_fit)
        # again subtract 1 to keep consistent with the fitting
        y, x = idx.I .- 1
        z_fit[idx] = a*x + b*y + c*x*y + d*x^2 + e*y^2 + f
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
